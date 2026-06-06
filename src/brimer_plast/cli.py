"""Brimer-PLAST command-line interface."""

from __future__ import annotations

import hashlib
import subprocess
from datetime import datetime
from pathlib import Path

import typer

from brimer_plast.genome import reverse_complement
from brimer_plast.log_config import configure_logging, get_logger
from brimer_plast.pdf_report import build_pdf_report
from brimer_plast.pipeline import PipelineResult, dump_debug_info, run_pipeline
from brimer_plast.primer import (
    PRIMER_NUM_RETURN,
    PRIMER_OPT_SIZE,
    PRIMER_MIN_SIZE,
    PRIMER_MAX_SIZE,
    PRIMER_OPT_TM,
    PRIMER_MIN_TM,
    PRIMER_MAX_TM,
    PRIMER_MIN_GC,
    PRIMER_MAX_GC,
    PRIMER_PRODUCT_MIN,
    PRIMER_PRODUCT_MAX,
    DEFAULT_PRIMER_ARGS,
)

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


def _run_for_target(
    *,
    target_key: str,
    target_type: str,  # "gene" or "transcript"
    genome: Path,
    annotations: Path,
    disable_junction_overlap: bool,
    num_return: int,
    min_tm: float,
    max_tm: float,
    opt_tm: float,
    min_size: int,
    max_size: int,
    opt_size: int,
    min_gc: float,
    max_gc: float,
    product_size_min: int,
    product_size_max: int,
    max_amplicon: int,
    tsv: bool,
    verbose: int,
    pdf_path: str | None,
) -> None:
    """Run the full Brimer-PLAST pipeline for a single target."""
    log = get_logger()

    target_gene = target_key if target_type == "gene" else None
    target_transcript = target_key if target_type == "transcript" else None

    # Build primer3 args from CLI parameters
    primer_args = {
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

    try:
        result = run_pipeline(
            genome=genome,
            annotations=annotations,
            target_key=target_key,
            target_type=target_type,
            disable_junction_overlap=disable_junction_overlap,
            primer_args=primer_args,
            max_amplicon=max_amplicon,
        )
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except (RuntimeError, FileNotFoundError) as e:
        typer.echo(f"  tnBLAST error: {e}", err=True)
        raise typer.Exit(code=1)

    chains = result.chains
    locus = result.locus
    filtered = result.filtered_pairs
    flat_pairs = result.all_candidates

    # ── Output results ─────────────────────────────────────────────────────
    typer.echo("")
    if filtered:
        if tsv:
            typer.echo(
                "pair_name\tforward_seq\treverse_seq\treverse_rc\tforward_tm\t"
                "reverse_tm\tforward_gc\treverse_gc\tproduct_size"
            )
            for pair in filtered:
                rc_rev = reverse_complement(pair.reverse_seq or "")
                typer.echo(
                    f"{pair.pair_name}\t{pair.forward_seq}\t{pair.reverse_seq}\t{rc_rev}\t"
                    f"{pair.forward_tm:.1f}\t{pair.reverse_tm:.1f}\t"
                    f"{pair.forward_gc:.0f}\t{pair.reverse_gc:.0f}\t"
                    f"{pair.product_size}"
                )
        else:
            typer.echo(
                f"{'Pair Name':<20} {'Forward':<28} {'Tm(°C)':<8} {'%GC':<5} "
                f"{'Reverse':<28} {'Tm(°C)':<8} {'%GC':<5} {'Size':<6}"
            )
            typer.echo("-" * 120)
            for pair in filtered:
                rc_rev = reverse_complement(pair.reverse_seq or "")
                typer.echo(
                    f"{pair.pair_name:<20} {pair.forward_seq:<28} {pair.forward_tm:<8.1f} "
                    f"{pair.forward_gc:<5.0f} {pair.reverse_seq:<28} "
                    f"{pair.reverse_tm:<8.1f} {pair.reverse_gc:<5.0f} "
                    f"{pair.product_size:<6}"
                )
                typer.echo(
                    f"{'':<20} {'':<28} {'':<8} {'':<5} {f'({rc_rev})':<28}"
                )
    else:
        typer.echo("No specificity-filtered primer pairs to display.")

    # ── Debug dump (independent of PDF generation) ────────────────────────
    if verbose >= 2:
        dump_debug_info(result)

    # ── Generate PDF report ────────────────────────────────────────────────
    if pdf_path is not None:
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
    target_gene: list[str] | None = typer.Option(
        None,
        "--target-gene",
        help="Target gene name (e.g. GAPDH). Repeat for multiple targets. "
        "One of --target-gene or --target-transcript is required.",
    ),
    target_transcript: list[str] | None = typer.Option(
        None,
        "--target-transcript",
        help="Target transcript ID (e.g. NM_001289746.1). Repeat for multiple "
        "targets. One of --target-gene or --target-transcript is required.",
    ),
    disable_junction_overlap: bool = typer.Option(
        False,
        "--disable-junction-overlap",
        help="Allow primers that do not span an exon-exon junction. "
        "Use this for genomic PCR rather than qRT-PCR.",
    ),
    num_return: int = typer.Option(
        PRIMER_NUM_RETURN,
        "--num-return",
        "-n",
        help="Number of candidate primer pairs to design (before filtering).",
    ),
    min_tm: float = typer.Option(
        PRIMER_MIN_TM,
        "--min-tm",
        help="Minimum primer melting temperature.",
    ),
    max_tm: float = typer.Option(
        PRIMER_MAX_TM,
        "--max-tm",
        help="Maximum primer melting temperature.",
    ),
    opt_tm: float = typer.Option(
        PRIMER_OPT_TM,
        "--opt-tm",
        help="Optimal primer melting temperature.",
    ),
    min_size: int = typer.Option(
        PRIMER_MIN_SIZE,
        "--min-size",
        help="Minimum primer length.",
    ),
    max_size: int = typer.Option(
        PRIMER_MAX_SIZE,
        "--max-size",
        help="Maximum primer length.",
    ),
    opt_size: int = typer.Option(
        PRIMER_OPT_SIZE,
        "--opt-size",
        help="Optimal primer length.",
    ),
    min_gc: float = typer.Option(
        PRIMER_MIN_GC,
        "--min-gc",
        help="Minimum primer GC content (percent).",
    ),
    max_gc: float = typer.Option(
        PRIMER_MAX_GC,
        "--max-gc",
        help="Maximum primer GC content (percent).",
    ),
    product_size_min: int = typer.Option(
        PRIMER_PRODUCT_MIN,
        "--product-min",
        help="Minimum PCR product size (bp).",
    ),
    product_size_max: int = typer.Option(
        PRIMER_PRODUCT_MAX,
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
    output_pdf: list[Path] | None = typer.Option(
        None,
        "--output-pdf",
        help="Write PDF report to this path (implies PDF generation). "
        "Repeat once per target, or omit for auto-generated names.",
        exists=False,
        dir_okay=False,
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

    # ── Validate PDF options ────────────────────────────────────────────────
    if output_pdf and no_pdf:
        typer.echo(
            "Error: --output-pdf and --no-pdf are mutually exclusive.",
            err=True,
        )
        raise typer.Exit(code=1)

    # ── Determine the list of targets ──────────────────────────────────────
    if target_gene:
        targets: list[tuple[str, str]] = [("gene", g) for g in (target_gene or [])]
    else:
        targets = [("transcript", t) for t in (target_transcript or [])]

    n_targets = len(targets)

    if output_pdf and len(output_pdf) != n_targets:
        typer.echo(
            f"Error: Number of --output-pdf values ({len(output_pdf)}) must match "
            f"the number of targets ({n_targets}).",
            err=True,
        )
        raise typer.Exit(code=1)

    genome_md5 = calculate_md5(genome)
    annotations_md5 = calculate_md5(annotations)
    version_str = get_git_version()

    # ── Loop over each target (independent invocation) ─────────────────
    for idx, (target_type, target_key) in enumerate(targets):
        typer.echo(
            f"\n{'=' * 60}",
        )
        typer.echo(
            f"  {target_type}: {target_key}  ({idx + 1} of {n_targets})",
        )
        typer.echo(
            f"{'=' * 60}",
        )

        # Determine PDF path for this target
        if no_pdf:
            pdf_path: str | None = None
        elif output_pdf:
            pdf_path = str(output_pdf[idx])
        else:
            slug = target_key.replace("|", "_").replace("/", "_")
            pdf_path = f"brimer_plast_{slug}_{datetime.now():%Y%m%d_%H%M%S}.pdf"

        _run_for_target(
            target_key=target_key,
            target_type=target_type,
            genome=genome,
            annotations=annotations,
            disable_junction_overlap=disable_junction_overlap,
            num_return=num_return,
            min_tm=min_tm,
            max_tm=max_tm,
            opt_tm=opt_tm,
            min_size=min_size,
            max_size=max_size,
            opt_size=opt_size,
            min_gc=min_gc,
            max_gc=max_gc,
            product_size_min=product_size_min,
            product_size_max=product_size_max,
            max_amplicon=max_amplicon,
            tsv=tsv,
            verbose=verbose,
            pdf_path=pdf_path,
        )
