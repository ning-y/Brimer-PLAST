"""Brimer-PLAST command-line interface."""

import hashlib
import logging
import os
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

import typer

from brimer_plast.filter import filter_specific_pairs, run_tnblast, write_assay_file
from brimer_plast.genome import (
    _exons_in_template_order,
    _reverse_complement,
    build_transcript_to_gene_map,
    build_transcriptome_fasta,
    genomic_range_to_fragments,
    get_gene_locus,
    get_target_information,
    template_to_genomic,
)
from brimer_plast.log_config import configure_logging, get_logger
from brimer_plast.models import ConservedExonChain, GeneLocus, GenomicFragment, PrimerPair
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
        result = subprocess.run(
            ["git", "describe", "--always", "--dirty"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "0.1.0"


def _dump_debug_info(
    log: logging.Logger,
    chains: list[ConservedExonChain],
    locus: GeneLocus | None,
    filtered_pairs: list[PrimerPair],
    all_flat_pairs: list[PrimerPair],
) -> None:
    """Log detailed debug information at DEBUG level."""
    log.debug("")
    log.debug("=" * 70)
    log.debug("DEBUG DUMP")
    log.debug("=" * 70)

    if locus is not None:
        log.debug("")
        log.debug("--- Locus ---")
        log.debug(f"  Gene: {locus.gene_name}")
        log.debug(f"  Seqid: {locus.seqid}")
        log.debug(f"  Strand: {locus.strand}")
        log.debug(f"  Genomic range: {locus.min_start:,} - {locus.max_end:,}")
        log.debug(f"  Transcripts ({len(locus.transcripts)}):")
        for tid, exons in sorted(locus.transcripts.items()):
            log.debug(f"    {tid}: {len(exons)} exons")
            for ex in exons:
                log.debug(f"      {ex.seqid}:{ex.start:,}-{ex.end:,} ({ex.strand})")

        g_min, g_max = locus.min_start, locus.max_end
        p_coords = []
        for p in filtered_pairs:
            for frag in (
                list(p.primer3_forward_fragments)
                + list(p.primer3_reverse_fragments)
                + list(p.tnblast_forward_fragments)
                + list(p.tnblast_reverse_fragments)
            ):
                p_coords.extend([frag.start, frag.end])
        if p_coords:
            pmin, pmax = min(p_coords), max(p_coords)
            pad = max(500, int((pmax - pmin) * 0.4))
            v_min = max(g_min, pmin - pad)
            v_max = min(g_max, pmax + pad)
        else:
            v_min, v_max = g_min, g_max
        log.debug(f"  Zoom range: {v_min:,} - {v_max:,}")
        log.debug(f"  Overview: {g_min:,} - {g_max:,}")

    log.debug("")
    log.debug("--- Conserved Exon Chains ---")
    for chain in chains:
        ordered = _exons_in_template_order(chain.exons)
        log.debug(f"  Chain: {chain.id}")
        log.debug(f"    Template length: {len(chain.template)} bp")
        log.debug(f"    Exons ({len(ordered)}):")
        cumulative = 1
        for ex in ordered:
            exon_len = ex.end - ex.start + 1
            if ex.strand == "-":
                coord_str = f"{ex.seqid}:{ex.end:,}-{ex.start:,} (-)"
            else:
                coord_str = f"{ex.seqid}:{ex.start:,}-{ex.end:,} (+)"
            log.debug(
                f"      {coord_str}  [{exon_len} bp,"
                f" template {cumulative}-{cumulative + exon_len - 1}]"
            )
            cumulative += exon_len
        log.debug(f"    Junction positions (1-based): {chain.junction_positions_1based}")
        log.debug(f"    Required junction positions: {chain.required_junction_positions_1based}")

    log.debug("")
    log.debug("--- Filtered Primer Pairs ---")
    for pair in filtered_pairs:
        pnum = pair.pair_number or "?"
        log.debug("")
        log.debug(f"  Pair {pnum} (chain: {pair.chain_id})")
        log.debug(f"    Product size: {pair.product_size}  Penalty: {pair.pair_penalty}")
        log.debug(f"    Forward ({pair.forward_len} bp): {pair.forward_seq}")
        log.debug(f"      Primer3 (blue): {[(f.seqid, f.start, f.end) for f in pair.primer3_forward_fragments]}")
        log.debug(f"      tnBLAST (red):  {[(f.seqid, f.start, f.end) for f in pair.tnblast_forward_fragments]}")
        log.debug(f"    Reverse ({pair.reverse_len} bp): {pair.reverse_seq}")
        log.debug(f"      Primer3 (blue): {[(f.seqid, f.start, f.end) for f in pair.primer3_reverse_fragments]}")
        log.debug(f"      tnBLAST (red):  {[(f.seqid, f.start, f.end) for f in pair.tnblast_reverse_fragments]}")
        if pair.forward_start is not None:
            log.debug(f"      Forward template 1-based: {pair.forward_start + 1}-{pair.forward_start + pair.forward_len}")
        if pair.reverse_start is not None:
            r1 = pair.reverse_start - pair.reverse_len + 2
            r2 = pair.reverse_start + 1
            log.debug(f"      Reverse template 1-based: {r1}-{r2}")

    log.debug("")
    log.debug("--- All Candidate Pairs (includes filtered-out) ---")
    for i, pair in enumerate(all_flat_pairs, start=1):
        mark = " [FILTERED IN]" if pair.pair_number is not None else ""
        log.debug(f"  pair_{i}: {pair.forward_seq} / {pair.reverse_seq}  chain={pair.chain_id}{mark}")
        log.debug(f"    P3 F: {pair.primer3_forward_fragments}")
        log.debug(f"    P3 R: {pair.primer3_reverse_fragments}")
        log.debug(f"    tn F: {pair.tnblast_forward_fragments}")
        log.debug(f"    tn R: {pair.tnblast_reverse_fragments}")

    log.debug("=" * 70)
    log.debug("")


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
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        count=True,
        help="Increase verbosity. -v for pipeline progress, -vv adds per-pair "
        "fragment-list details and template coordinates.",
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
    configure_logging(verbose)
    log = get_logger()

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
    log.info("Reading genome and annotations...")
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

    log.info(
        f"  Found {len(chains)} conserved exon chain(s)"
    )
    for chain in chains:
        log.info(
            f"    {chain.id}: {len(chain.template)} bp template, "
            f"{len(chain.junction_positions_1based)} junction(s)"
        )

    # ── Step 2: Design candidate primers for each chain ────────────────────
    log.info("Designing primers with primer3...")
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
            log.info(
                f"    {chain.id}: {len(candidate_pairs)} candidate pair(s)"
            )
            all_candidates.append((chain, candidate_pairs))
        else:
            log.warning(
                f"    {chain.id}: no candidate primers could be designed."
            )

    if not all_candidates:
        typer.echo("  No candidate primers could be designed for any chain.", err=True)
        raise typer.Exit(code=1)

    # ── Step 3: Flat list for tnBLAST (merge all chains) ───────────────────
    flat_pairs = []
    for chain, pairs in all_candidates:
        flat_pairs.extend(pairs)
    log.info(
        f"  Total: {len(flat_pairs)} candidate pair(s) across all chains."
    )

    # ── Step 4: Filter with tnBLAST ────────────────────────────────────────
    log.info("Running tnBLAST specificity filter...")

    with tempfile.TemporaryDirectory() as tmp_dir:
        assay_path = os.path.join(tmp_dir, "assays.txt")
        write_assay_file(flat_pairs, assay_path)

        # Build transcriptome from annotation
        typer.echo("  Building transcriptome from annotations...", err=True)
        transcriptome_path = os.path.join(tmp_dir, "transcriptome.fa")
        transcript_exon_map = build_transcriptome_fasta(
            genome, annotations, transcriptome_path
        )

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
            log.info("  tnBLAST genome scan complete")
            run_tnblast(
                assay_path,
                transcriptome_path,
                max_amplicon=max_amplicon,
                min_tm=min_tm,
                max_tm=max_tm,
                output_path=transcriptome_out,
            )
            log.info("  tnBLAST transcriptome scan complete")
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

        # Build chain map for exon lookup
        chain_map = {c.id: c for c in chains}

        # Compute genomic fragment lists for all pairs
        for i, pair in enumerate(flat_pairs, start=1):
            name = f"pair_{i}"

            # ── Primer3-derived fragments (Panel A, blue) ──
            if pair.forward_start is not None and pair.chain_id in chain_map:
                exons = chain_map[pair.chain_id].exons
                pair.primer3_forward_fragments = template_to_genomic(
                    pair.forward_start, pair.forward_len or 20, exons
                )
                # PRIMER_RIGHT_n returns (3' end, length) — 0-based.
                # Convert to 5' start for template_to_genomic.
                rev_5prime = pair.reverse_start - pair.reverse_len + 1
                pair.primer3_reverse_fragments = template_to_genomic(
                    rev_5prime, pair.reverse_len or 20, exons
                )

            # ── tnBLAST-derived fragments (Panel B, red) ──
            if pair.chain_id in chain_map:
                exons = chain_map[pair.chain_id].exons
                f_len = pair.forward_len or 20
                r_len = pair.reverse_len or 20

                if name in genome_amplicons and genome_amplicons[name]:
                    hit = genome_amplicons[name][0]
                    if locus and locus.strand == "-":
                        # tnBLAST amplicon range is 0-based inclusive.
                        # For a negative-strand gene, the forward primer
                        # sits on the minus strand; its reverse complement
                        # appears on the plus strand at the amplicon end.
                        # Convert to 1-based inclusive for
                        # genomic_range_to_fragments.
                        f_g_start = hit.amplicon_end - f_len + 2
                        f_g_end = hit.amplicon_end + 1
                        r_g_start = hit.amplicon_start + 1
                        r_g_end = hit.amplicon_start + r_len
                    else:
                        f_g_start = hit.amplicon_start + 1
                        f_g_end = hit.amplicon_start + f_len
                        r_g_start = hit.amplicon_end - r_len + 2
                        r_g_end = hit.amplicon_end + 1

                    pair.tnblast_forward_fragments = genomic_range_to_fragments(
                        f_g_start, f_g_end, exons
                    )
                    pair.tnblast_reverse_fragments = genomic_range_to_fragments(
                        r_g_start, r_g_end, exons
                    )

                elif name in transcriptome_amplicons and transcriptome_amplicons[name]:
                    hit = transcriptome_amplicons[name][0]
                    # tnBLAST transcriptome coordinates are 0-based inclusive
                    # in the hit transcript's coordinate space.  Map through
                    # that specific transcript's exon set so the red marker
                    # shows where tnBLAST independently found the primer.
                    tr_exons = transcript_exon_map.get(hit.seqid)
                    if tr_exons is not None:
                        pair.tnblast_forward_fragments = template_to_genomic(
                            hit.amplicon_start, f_len, tr_exons
                        )
                        pair.tnblast_reverse_fragments = template_to_genomic(
                            hit.amplicon_end - r_len + 1, r_len, tr_exons
                        )

    if not filtered:
        log.warning("  No primer pairs passed tnBLAST specificity filter.")
    else:
        log.info(
            f"  {len(filtered)}/{len(flat_pairs)} pairs passed filtering."
        )

    # Assign pair numbers
    for i, pair in enumerate(filtered, start=1):
        pair.pair_number = i

    # ── Step 5: Output results ─────────────────────────────────────────────
    typer.echo("")
    if filtered:
        if tsv:
            typer.echo(
                "pair\tforward_seq\treverse_seq\treverse_rc\tforward_tm\t"
                "reverse_tm\tforward_gc\treverse_gc\tproduct_size"
            )
            for i, pair in enumerate(filtered, start=1):
                rc_rev = _reverse_complement(pair.reverse_seq or "")
                typer.echo(
                    f"{i}\t{pair.forward_seq}\t{pair.reverse_seq}\t{rc_rev}\t"
                    f"{pair.forward_tm:.1f}\t{pair.reverse_tm:.1f}\t"
                    f"{pair.forward_gc:.0f}\t{pair.reverse_gc:.0f}\t"
                    f"{pair.product_size}"
                )
        else:
            typer.echo(
                f"{'Pair':<6} {'Forward':<28} {'Tm(°C)':<8} {'%GC':<5} "
                f"{'Reverse':<28} {'Tm(°C)':<8} {'%GC':<5} {'Size':<6}"
            )
            typer.echo("-" * 100)
            for i, pair in enumerate(filtered, start=1):
                rc_rev = _reverse_complement(pair.reverse_seq or "")
                typer.echo(
                    f"{i:<6} {pair.forward_seq:<28} {pair.forward_tm:<8.1f} "
                    f"{pair.forward_gc:<5.0f} {pair.reverse_seq:<28} "
                    f"{pair.reverse_tm:<8.1f} {pair.reverse_gc:<5.0f} "
                    f"{pair.product_size:<6}"
                )
                typer.echo(
                    f"{'':<6} {'':<28} {'':<8} {'':<5} {f'({rc_rev})':<28}"
                )
    else:
        typer.echo("No specificity-filtered primer pairs to display.")

    # ── Debug dump (independent of PDF generation) ────────────────────────
    if verbose >= 2:
        _dump_debug_info(log, chains, locus, filtered, flat_pairs)

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
