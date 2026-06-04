"""Brimer-PLAST core pipeline.

Extracts the primer design and specificity-filtering logic from the CLI
into a reusable function that returns a structured result.
"""

from __future__ import annotations

import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from brimer_plast.filter import filter_specific_pairs
from brimer_plast.tnblast import (
    _parse_tnblast_amplicons,
    run_tnblast,
    write_assay_file,
)
from brimer_plast.genome import (
    build_transcript_to_gene_map,
    build_transcriptome_fasta,
    exons_in_template_order,
    genomic_range_to_fragments,
    get_gene_locus,
    get_target_information,
    template_to_genomic,
)
from brimer_plast.models import (
    ConservedExonChain,
    ExonInfo,
    GeneLocus,
    PrimerPair,
)
from brimer_plast.primer import design_primers


@dataclass
class PipelineResult:
    """Structured result from a single pipeline run."""

    chains: list[ConservedExonChain] = field(default_factory=list)
    locus: GeneLocus | None = None
    filtered_pairs: list[PrimerPair] = field(default_factory=list)
    all_candidates: list[PrimerPair] = field(default_factory=list)

    @property
    def has_results(self) -> bool:
        """Whether any filtered pairs survived tnBLAST screening."""
        return len(self.filtered_pairs) > 0


def run_pipeline(
    *,
    genome: Path,
    annotations: Path,
    target_key: str,
    target_type: str,  # "gene" or "transcript"
    disable_junction_overlap: bool,
    primer_args: dict[str, Any],
    max_amplicon: int = 2000,
) -> PipelineResult:
    """Run the full Brimer-PLAST design + filter pipeline for one target.

    Handles genome/annotation parsing, primer3 design, tnBLAST I/O, and
    specificity filtering.  The caller is responsible for CLI output and
    PDF generation.

    Args:
        genome: Path to genome FASTA.
        annotations: Path to annotation GTF.
        target_key: The gene name or transcript ID.
        target_type: ``"gene"`` or ``"transcript"``.
        disable_junction_overlap: If True, do not enforce junction-spanning.
        primer_args: primer3 global args (merged over defaults).
        max_amplicon: Maximum tnBLAST amplicon search length.

    Returns:
        A :class:`PipelineResult` with chains, filtered pairs, etc.

    Raises:
        ValueError: if the target is not found or has no usable data.
        RuntimeError: if tnBLAST fails.
        FileNotFoundError: if tnBLAST is not on PATH.
    """
    log = logging.getLogger("brimer_plast.pipeline")

    target_gene = target_key if target_type == "gene" else None
    target_transcript = target_key if target_type == "transcript" else None

    # ── Step 1: Extract conserved exon chains ──────────────────────────────
    chains = get_target_information(
        fasta_path=genome,
        gtf_path=annotations,
        target_gene=target_gene,
        target_transcript=target_transcript,
    )
    locus = get_gene_locus(
        gtf_path=annotations,
        target_gene=target_gene,
        target_transcript=target_transcript,
    )

    if not disable_junction_overlap:
        chains_with_junctions = [c for c in chains if c.junction_positions_1based]
        if not chains_with_junctions:
            raise ValueError(
                "No multi-exon targets found and --disable-junction-overlap "
                "was not set.  The target gene may consist only of single-exon "
                "transcripts (e.g. mitochondrial genes).  Add "
                "--disable-junction-overlap to design primers without the "
                "junction-spanning requirement."
            )
        chains = chains_with_junctions

    log.info("Found %d conserved exon chain(s)", len(chains))
    for chain in chains:
        log.info(
            "  %s: %d bp template, %d junction(s)",
            chain.id, len(chain.template), len(chain.junction_positions_1based),
        )

    # ── Step 2: Design candidate primers for each chain ────────────────────
    log.info("Designing primers with primer3...")
    all_candidates: list[tuple[ConservedExonChain, list[PrimerPair]]] = []
    for chain in chains:
        if disable_junction_overlap:
            junction_positions: list[int] = []
            required_junction_positions: list[int] | None = None
        else:
            junction_positions = chain.junction_positions_1based
            required_junction_positions = (
                chain.required_junction_positions_1based or junction_positions
            )
        candidate_pairs = design_primers(
            chain.template,
            sequence_id=chain.id,
            chain_id=chain.id,
            global_args=primer_args,
            junction_positions=junction_positions,
            required_junction_positions=required_junction_positions,
        )
        if candidate_pairs:
            log.info("    %s: %d candidate pair(s)", chain.id, len(candidate_pairs))
            all_candidates.append((chain, candidate_pairs))
        else:
            log.warning("    %s: no candidate primers could be designed.", chain.id)

    if not all_candidates:
        raise ValueError("No candidate primers could be designed for any chain.")

    flat_pairs = [p for _, pairs in all_candidates for p in pairs]
    log.info("Total: %d candidate pair(s) across all chains.", len(flat_pairs))

    # ── Step 3: Filter with tnBLAST ────────────────────────────────────────
    log.info("Running tnBLAST specificity filter...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        assay_path = os.path.join(tmp_dir, "assays.txt")
        write_assay_file(flat_pairs, assay_path)

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
                min_tm=primer_args.get("PRIMER_MIN_TM", 57.0),
                max_tm=primer_args.get("PRIMER_MAX_TM", 63.0),
                output_path=genome_out,
            )
            transcriptome_out = os.path.join(tmp_dir, "tntblast_transcriptome.txt")
            log.info("  tnBLAST genome scan complete")
            run_tnblast(
                assay_path,
                transcriptome_path,
                max_amplicon=max_amplicon,
                min_tm=primer_args.get("PRIMER_MIN_TM", 57.0),
                max_tm=primer_args.get("PRIMER_MAX_TM", 63.0),
                output_path=transcriptome_out,
            )
            log.info("  tnBLAST transcriptome scan complete")
        except (RuntimeError, FileNotFoundError):
            raise

        genome_amplicons = _parse_tnblast_amplicons(genome_out)
        transcriptome_amplicons = _parse_tnblast_amplicons(transcriptome_out)

        t2g = build_transcript_to_gene_map(annotations)
        transcriptome_targets: dict[str, list[str]] = {}
        for name, hits in transcriptome_amplicons.items():
            genes = set()
            for hit in hits:
                tid = hit.seqid
                gene = t2g.get(tid, tid)
                genes.add(gene)
            transcriptome_targets[name] = list(genes)

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

        # Assign pair numbers
        for i, pair in enumerate(filtered, start=1):
            pair.pair_number = i

        # ── Step 4: Compute genomic fragment lists ─────────────────────────
        chain_map = {c.id: c for c in chains}

        for i, pair in enumerate(flat_pairs, start=1):
            name = f"pair_{i}"

            if (
                pair.forward_start is not None
                and pair.reverse_start is not None
                and pair.chain_id in chain_map
            ):
                exons = chain_map[pair.chain_id].exons
                pair.primer3_forward_fragments = template_to_genomic(
                    pair.forward_start, pair.forward_len or 20, exons
                )
                rev_5prime = pair.reverse_start - (pair.reverse_len or 20) + 1
                pair.primer3_reverse_fragments = template_to_genomic(
                    rev_5prime, pair.reverse_len or 20, exons
                )

            if pair.chain_id in chain_map:
                exons = chain_map[pair.chain_id].exons
                f_len = pair.forward_len or 20
                r_len = pair.reverse_len or 20

                if name in genome_amplicons and genome_amplicons[name]:
                    hit = genome_amplicons[name][0]
                    if locus and locus.strand == "-":
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
                    tr_exons = transcript_exon_map.get(hit.seqid)
                    if tr_exons is not None:
                        pair.tnblast_forward_fragments = template_to_genomic(
                            hit.amplicon_start, f_len, tr_exons
                        )
                        pair.tnblast_reverse_fragments = template_to_genomic(
                            hit.amplicon_end - r_len + 1, r_len, tr_exons
                        )

    if not filtered:
        log.warning("No primer pairs passed tnBLAST specificity filter.")
    else:
        log.info("%d/%d pairs passed filtering.", len(filtered), len(flat_pairs))

    return PipelineResult(
        chains=chains,
        locus=locus,
        filtered_pairs=filtered,
        all_candidates=flat_pairs,
    )


def dump_debug_info(
    result: PipelineResult,
    log: logging.Logger | None = None,
) -> None:
    """Log detailed debug information about a pipeline result at DEBUG level.

    Args:
        result: The pipeline result to dump.
        log: Logger to use.  If None, uses ``brimer_plast.pipeline``.
    """
    if log is None:
        log = logging.getLogger("brimer_plast.pipeline")

    chains = result.chains
    locus = result.locus
    filtered_pairs = result.filtered_pairs
    all_flat_pairs = result.all_candidates

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
        ordered = exons_in_template_order(chain.exons)
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
        if pair.forward_start is not None and pair.forward_len is not None:
            log.debug(f"      Forward template 1-based: {pair.forward_start + 1}-{pair.forward_start + pair.forward_len}")
        if pair.reverse_start is not None and pair.reverse_len is not None:
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
