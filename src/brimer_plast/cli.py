"""Brimer-PLAST command-line interface."""

import hashlib
import os
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

import typer

from brimer_plast.filter import filter_specific_pairs, run_tnblast, write_assay_file
from brimer_plast.genome import (
    build_transcript_to_gene_map,
    build_transcriptome_fasta,
    get_gene_locus,
    get_target_information,
    template_to_genomic,
)
from brimer_plast.models import ConservedExonChain, GeneLocus, PrimerPair
from brimer_plast.pdf_report import build_pdf_report
from brimer_plast.primer import DEFAULT_PRIMER_ARGS, design_primers

app = typer.Typer(
    name="brimer-plast",
    help="Design qRT-PCR primers (exon-exon junction-spanning) from a genome "
    "and annotation using primer3, then filter for specificity with tnBLAST.",
)


def calculate_md5(path: Path) -> str:
    """Calculate MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_git_version() -> str:
    """Get the current git version string if available."""
    try:
        # Check if we're in a git repo
        result = subprocess.run(
            ["git", "describe", "--always", "--dirty"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "0.1.0"


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    genome: Path = typer.Option(
        ...,
        "--genome",
        "-g",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Genome FASTA file.",
    ),
    annotations: Path = typer.Option(
        ...,
        "--annotations",
        "-a",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Gene annotation GTF file.",
    ),
    target_gene: Optional[str] = typer.Option(
        None,
        "--target-gene",
        help="Target gene name (e.g. GAPDH). One of --target-gene or "
        "--target-transcript is required.",
    ),
    target_transcript: Optional[str] = typer.Option(
        None,
        "--target-transcript",
        help="Target transcript ID (e.g. NM_001289746.1). One of --target-gene "
        "or --target-transcript is required.",
    ),
    disable_junction_overlap: bool = typer.Option(
        False,
        "--disable-junction-overlap",
        help="Allow primers that do not span an exon-exon junction. "
        "Use this for genomic PCR rather than qRT-PCR.",
    ),
    num_return: int = typer.Option(
        10,
        "--num-return",
        "-n",
        help="Number of candidate primer pairs to design (before filtering).",
    ),
    min_tm: float = typer.Option(
        57.0,
        "--min-tm",
        help="Minimum primer melting temperature.",
    ),
    max_tm: float = typer.Option(
        63.0,
        "--max-tm",
        help="Maximum primer melting temperature.",
    ),
    opt_tm: float = typer.Option(
        60.0,
        "--opt-tm",
        help="Optimal primer melting temperature.",
    ),
    min_size: int = typer.Option(
        18,
        "--min-size",
        help="Minimum primer length.",
    ),
    max_size: int = typer.Option(
        25,
        "--max-size",
        help="Maximum primer length.",
    ),
    opt_size: int = typer.Option(
        20,
        "--opt-size",
        help="Optimal primer length.",
    ),
    min_gc: float = typer.Option(
        40.0,
        "--min-gc",
        help="Minimum primer GC content (percent).",
    ),
    max_gc: float = typer.Option(
        60.0,
        "--max-gc",
        help="Maximum primer GC content (percent).",
    ),
    product_size_min: int = typer.Option(
        100,
        "--product-min",
        help="Minimum PCR product size (bp).",
    ),
    product_size_max: int = typer.Option(
        400,
        "--product-max",
        help="Maximum PCR product size (bp).",
    ),
    max_amplicon: int = typer.Option(
        2000,
        "--max-amplicon",
        help="Maximum tnBLAST amplicon search length.",
    ),
    tsv: bool = typer.Option(
        False,
        "--tsv",
        help="Output results as tab-separated values (machine-readable).",
    ),
    no_pdf: bool = typer.Option(
        False,
        "--no-pdf",
        help="Suppress PDF report generation.",
    ),
) -> None:
    """Design primers for a target and filter for specificity.

    By default, at least one primer in each pair must span an exon-exon
    junction (qRT-PCR mode).  Use --disable-junction-overlap for
    genomic PCR.
    """
    # ── Validate target arguments ──────────────────────────────────────────
    if not target_gene and not target_transcript:
        typer.echo(
            "Error: Provide either --target-gene or --target-transcript.",
            err=True,
        )
        raise typer.Exit(code=1)
    if target_gene and target_transcript:
        typer.echo(
            "Error: Provide --target-gene or --target-transcript, not both.",
            err=True,
        )
        raise typer.Exit(code=1)

    # ── Step 1: Extract conserved exon chains ──────────────────────────────
    typer.echo("Reading genome and annotations...", err=True)
    try:
        chains: list[ConservedExonChain] = get_target_information(
            fasta_path=genome,
            gtf_path=annotations,
            target_gene=target_gene,
            target_transcript=target_transcript,
        )
        locus: Optional[GeneLocus] = get_gene_locus(
            gtf_path=annotations,
            target_gene=target_gene,
            target_transcript=target_transcript,
        )
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)

    if not disable_junction_overlap:
        # Drop junctionless chains (single-exon transcripts) — they can't
        # satisfy the junction-spanning requirement.
        chains_with_junctions = [c for c in chains if c.junction_positions_1based]
        if not chains_with_junctions:
            typer.echo(
                "  Error: No multi-exon targets found and --disable-junction-overlap "
                "was not set.  The target gene may consist only of single-exon "
                "transcripts (e.g. mitochondrial genes).  Add "
                "--disable-junction-overlap to design primers without the "
                "junction-spanning requirement.",
                err=True,
            )
            raise typer.Exit(code=1)
        chains = chains_with_junctions

    typer.echo(
        f"  Found {len(chains)} conserved exon chain(s)", err=True
    )
    for chain in chains:
        typer.echo(
            f"    {chain.id}: {len(chain.template)} bp template, "
            f"{len(chain.junction_positions_1based)} junction(s)",
            err=True,
        )

    # ── Step 2: Design candidate primers for each chain ────────────────────
    typer.echo("Designing primers with primer3...", err=True)
    global_args = {
        **DEFAULT_PRIMER_ARGS,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_PRODUCT_SIZE_RANGE": f"{product_size_min}-{product_size_max}",
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": min_tm,
        "PRIMER_MAX_TM": max_tm,
        "PRIMER_OPT_SIZE": opt_size,
        "PRIMER_MIN_SIZE": min_size,
        "PRIMER_MAX_SIZE": max_size,
        "PRIMER_MIN_GC": min_gc,
        "PRIMER_MAX_GC": max_gc,
    }

    all_candidates: list[tuple[ConservedExonChain, list]] = []
    for chain in chains:
        if disable_junction_overlap:
            junction_positions = []
            required_junction_positions = None
        else:
            junction_positions = chain.junction_positions_1based
            # Use required_junction_positions_1based if set (unique junctions),
            # otherwise default to all junctions for post-filtering
            required_junction_positions = (
                chain.required_junction_positions_1based or junction_positions
            )
        candidate_pairs = design_primers(
            chain.template,
            sequence_id=chain.id,
            chain_id=chain.id,
            global_args=global_args,
            junction_positions=junction_positions,
            required_junction_positions=required_junction_positions,
        )
        if candidate_pairs:
            typer.echo(
                f"    {chain.id}: {len(candidate_pairs)} candidate pair(s)",
                err=True,
            )
            all_candidates.append((chain, candidate_pairs))
        else:
            typer.echo(
                f"    {chain.id}: no candidate primers could be designed.",
                err=True,
            )

    if not all_candidates:
        typer.echo("  No candidate primers could be designed for any chain.", err=True)
        raise typer.Exit(code=1)

    # ── Step 3: Flat list for tnBLAST (merge all chains) ───────────────────
    flat_pairs = []
    for chain, pairs in all_candidates:
        flat_pairs.extend(pairs)
    typer.echo(
        f"  Total: {len(flat_pairs)} candidate pair(s) across all chains.",
        err=True,
    )

    # ── Step 4: Filter with tnBLAST ────────────────────────────────────────
    typer.echo("Running tnBLAST specificity filter...", err=True)

    with tempfile.TemporaryDirectory() as tmp_dir:
        assay_path = os.path.join(tmp_dir, "assays.txt")
        write_assay_file(flat_pairs, assay_path)

        # Build transcriptome from annotation
        typer.echo("  Building transcriptome from annotations...", err=True)
        transcriptome_path = os.path.join(tmp_dir, "transcriptome.fa")
        build_transcriptome_fasta(genome, annotations, transcriptome_path)

        try:
            genome_out = os.path.join(tmp_dir, "tntblast_genome.txt")
            genome_counts = run_tnblast(
                assay_path,
                genome,
                max_amplicon=max_amplicon,
                min_tm=min_tm,
                max_tm=max_tm,
                output_path=genome_out,
            )
            transcriptome_out = os.path.join(tmp_dir, "tntblast_transcriptome.txt")
            run_tnblast(
                assay_path,
                transcriptome_path,
                max_amplicon=max_amplicon,
                min_tm=min_tm,
                max_tm=max_tm,
                output_path=transcriptome_out,
            )
        except (RuntimeError, FileNotFoundError) as e:
            typer.echo(f"  tnBLAST error: {e}", err=True)
            raise typer.Exit(code=1)

        # Re-parse genome tnBLAST output for Panel B amplicon coordinates
        from brimer_plast.filter import _parse_tnblast_amplicons

        genome_amplicons = _parse_tnblast_amplicons(genome_out)

        # Re-parse transcriptome tnBLAST output for mapping
        transcriptome_amplicons = _parse_tnblast_amplicons(transcriptome_out)

        # Map transcript IDs to gene names
        t2g = build_transcript_to_gene_map(annotations)
        transcriptome_targets: dict[str, list[str]] = {}
        for name, hits in transcriptome_amplicons.items():
            genes = set()
            for hit in hits:
                tid = hit.seqid
                gene = t2g.get(tid, tid)
                genes.add(gene)
            transcriptome_targets[name] = list(genes)

        resolved_gene: str
        if target_gene:
            resolved_gene = target_gene
        elif target_transcript:
            resolved_gene = t2g.get(target_transcript, target_transcript)
        else:
            resolved_gene = ""

        filtered = filter_specific_pairs(
            flat_pairs,
            genome_counts,
            transcriptome_targets,
            target_gene=resolved_gene,
            junction_mode=not disable_junction_overlap,
        )

        # Attach coordinates for Panel B (Sanity Check)
        for i, pair in enumerate(flat_pairs, start=1):
            name = f"pair_{i}"

            # Priority 1: Use direct genome hits (only in --disable-junction-overlap mode)
            if name in genome_amplicons and genome_amplicons[name]:
                hit = genome_amplicons[name][0]
                pair.tntblast_seqid = hit.seqid
                pair.tntblast_amplicon_start = hit.amplicon_start
                pair.tntblast_amplicon_end = hit.amplicon_end
            # Priority 2: Use mRNA hits mapped back to genome (works for junction mode)
            elif name in transcriptome_amplicons and transcriptome_amplicons[name]:
                hit = transcriptome_amplicons[name][0]  # first mRNA hit
                tid = hit.seqid
                if locus and tid in locus.transcripts:
                    exons = locus.transcripts[tid]
                    # Map mRNA range (1-based) back to genomic fragments
                    frags_start = template_to_genomic(hit.amplicon_start - 1, 1, exons)
                    frags_end = template_to_genomic(hit.amplicon_end - 1, 1, exons)
                    if frags_start and frags_end:
                        # Extract single genomic point from first fragment
                        f_start = frags_start[0]
                        f_end = frags_end[0]
                        pair.tntblast_seqid = f_start[0]
                        pair.tntblast_amplicon_start = min(f_start[1], f_end[1])
                        pair.tntblast_amplicon_end = max(f_start[2], f_end[2])
    if not filtered:
        typer.echo("  No primer pairs passed tnBLAST specificity filter.", err=True)
    else:
        typer.echo(
            f"  {len(filtered)}/{len(flat_pairs)} pairs passed filtering.",
            err=True,
        )

    # Assign pair numbers
    for i, pair in enumerate(filtered, start=1):
        pair.pair_number = i

    # ── Step 5: Output results ─────────────────────────────────────────────
    typer.echo("")
    if filtered:
        if tsv:
            typer.echo(
                "pair\tforward_seq\treverse_seq\tforward_tm\t"
                "reverse_tm\tforward_gc\treverse_gc\tproduct_size"
            )
            for i, pair in enumerate(filtered, start=1):
                typer.echo(
                    f"{i}\t{pair.forward_seq}\t{pair.reverse_seq}\t"
                    f"{pair.forward_tm:.1f}\t{pair.reverse_tm:.1f}\t"
                    f"{pair.forward_gc:.0f}\t{pair.reverse_gc:.0f}\t"
                    f"{pair.product_size}"
                )
        else:
            typer.echo(
                f"{'Pair':<6} {'Forward (5→3)':<28} {'Tm(°C)':<8} {'%GC':<5} "
                f"{'Reverse (5→3)':<28} {'Tm(°C)':<8} {'%GC':<5} {'Size':<6}"
            )
            typer.echo("-" * 100)
            for i, pair in enumerate(filtered, start=1):
                typer.echo(
                    f"{i:<6} {pair.forward_seq:<28} {pair.forward_tm:<8.1f} "
                    f"{pair.forward_gc:<5.0f} {pair.reverse_seq:<28} "
                    f"{pair.reverse_tm:<8.1f} {pair.reverse_gc:<5.0f} "
                    f"{pair.product_size:<6}"
                )
    else:
        typer.echo("No specificity-filtered primer pairs to display.")

    # ── Step 6: Generate PDF report ───────────────────────────────────────
    if not no_pdf:
        gene_slug = (target_gene or target_transcript or "unknown").replace("|", "_").replace("/", "_")
        pdf_path = f"brimer_plast_{gene_slug}_{datetime.now():%Y%m%d_%H%M%S}.pdf"
        typer.echo(f"\nGenerating PDF report: {pdf_path}", err=True)
        try:
            genome_md5 = calculate_md5(genome)
            annotations_md5 = calculate_md5(annotations)
            version_str = get_git_version()

            build_pdf_report(
                output_path=pdf_path,
                chains=chains,
                locus=locus,
                filtered_pairs=filtered,
                target_gene=target_gene,
                target_transcript=target_transcript,
                genome_path=str(genome),
                annotations_path=str(annotations),
                genome_md5=genome_md5,
                annotations_md5=annotations_md5,
                version_str=version_str,
                cli_args={
                    "target_gene": target_gene,
                    "target_transcript": target_transcript,
                    "disable_junction_overlap": disable_junction_overlap,
                    "num_return": num_return,
                    "min_tm": min_tm,
                    "max_tm": max_tm,
                    "opt_tm": opt_tm,
                    "min_size": min_size,
                    "max_size": max_size,
                    "opt_size": opt_size,
                    "min_gc": min_gc,
                    "max_gc": max_gc,
                    "product_min": product_size_min,
                    "product_max": product_size_max,
                    "max_amplicon": max_amplicon,
                },
            )
            typer.echo(f"  PDF written to {pdf_path}", err=True)
        except Exception as e:
            typer.echo(f"  PDF generation failed: {e}", err=True)
