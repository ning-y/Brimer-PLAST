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
from brimer_plast.tnblast import AmpliconHit, _parse_tnblast_amplicons, run_tnblast, write_assay_file
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
from brimer_plast.primer import design_primers, design_primers_dual_mode, dedup_and_prioritize


@dataclass
class PipelineResult:
    """Structured result from a single pipeline run."""

    chains: list[ConservedExonChain] = field(default_factory=list)
    locus: GeneLocus | None = None
    filtered_pairs: list[PrimerPair] = field(default_factory=list)
    all_candidates: list[PrimerPair] = field(default_factory=list)
    junction_candidates: list[PrimerPair] = field(default_factory=list)
    intron_candidates: list[PrimerPair] = field(default_factory=list)

    # Warnings accumulated during the run (e.g. fallback to per-transcript
    # chains, chain-specific primer design failures).  The CLI reprints
    # these at the end of a multi-target loop so they aren't lost in
    # scrolling output.
    warnings: list[str] = field(default_factory=list)

    @property
    def has_results(self) -> bool:
        """Whether any filtered pairs survived tnBLAST screening."""
        return len(self.filtered_pairs) > 0


# ── helper: tnBLAST fragment builder ───────────────────────────────────────


def _compute_tnblast_fragments(
    pair: PrimerPair,
    name: str,
    exons: list[ExonInfo],
    locus: GeneLocus | None,
    genome_amplicons: dict[str, list[AmpliconHit]],
    transcriptome_amplicons: dict[str, list[AmpliconHit]],
    transcript_exon_map: dict[str, list[ExonInfo]],
) -> None:
    """Populate tnBLAST-derived fragment lists on a PrimerPair.

    Uses the first hit from genome amplicons if available (preferred),
    otherwise falls back to transcriptome amplicons.
    """
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


def run_pipeline(
    *,
    genome: Path,
    annotations: Path,
    target_key: str,
    target_type: str,  # "gene" or "transcript"
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

    # Separate single-exon chains (no junctions, no introns possible)
    # from multi-exon chains that can use both design modes.
    pipeline_warnings: list[str] = []
    fallback_chains = [c for c in chains if c.fallback]
    if fallback_chains:
        fallback_ids = [c.id for c in fallback_chains]
        msg = (
            f"No conserved exon-exon junctions were found across transcripts of "
            f"{target_key}. Primer design fell back to per-transcript chains "
            f"({', '.join(fallback_ids)}). Primers from these chains may not "
            f"target all transcripts of the gene — verify manually, or use "
            f"--target-transcript for single-transcript design."
        )
        log.warning(msg)
        pipeline_warnings.append(msg)

    multi_exon_chains = [c for c in chains if len(c.exons) >= 2]
    single_exon_chains = [c for c in chains if len(c.exons) < 2]

    if single_exon_chains:
        single_ids = [c.id for c in single_exon_chains]
        msg = (
            f"Single-exon chain(s) {', '.join(single_ids)} have no exon-exon "
            f"junctions and no introns. No primers can be designed for these "
            f"chains."
        )
        log.warning(msg)
        pipeline_warnings.append(msg)

    log.info("Found %d conserved exon chain(s) (%d multi-exon, %d single-exon)",
              len(chains), len(multi_exon_chains), len(single_exon_chains))

    # ── Step 2: Design candidate primers (dual-mode) ──────────────────────
    log.info("Designing primers with primer3 (junction + intron modes)...")
    num_return = primer_args.get("PRIMER_NUM_RETURN", 50)
    junction_candidates: list[PrimerPair] = []
    intron_candidates: list[PrimerPair] = []

    for chain in multi_exon_chains:
        jpos = chain.junction_positions_1based
        rpos = chain.required_junction_positions_1based or jpos

        j_pairs, i_pairs = design_primers_dual_mode(
            chain.template,
            chain.exons,
            sequence_id=chain.id,
            chain_id=chain.id,
            global_args=primer_args,
            num_return=num_return,
            junction_positions=jpos,
            required_junction_positions=rpos,
        )
        if j_pairs:
            log.info("    %s: %d junction-mode candidate(s)", chain.id, len(j_pairs))
            junction_candidates.extend(j_pairs)
        if i_pairs:
            log.info("    %s: %d intron-mode candidate(s)", chain.id, len(i_pairs))
            intron_candidates.extend(i_pairs)
        if not j_pairs and not i_pairs:
            log.warning("    %s: no candidate primers could be designed.", chain.id)

    all_flat_pairs = junction_candidates + intron_candidates
    if not all_flat_pairs:
        msg = (
            f"No candidate primers could be designed for any chain. "
            f"The target may consist only of single-exon transcripts, "
            f"or the design constraints may be too strict."
        )
        log.warning(msg)
        pipeline_warnings.append(msg)
        return PipelineResult(
            chains=chains,
            locus=locus,
            filtered_pairs=[],
            all_candidates=[],
            junction_candidates=[],
            intron_candidates=[],
            warnings=pipeline_warnings,
        )

    log.info("Total: %d candidate pair(s) (%d junction, %d intron).",
              len(all_flat_pairs), len(junction_candidates), len(intron_candidates))

    # ── Step 2.5: Compute pair names ───────────────────────────────────────
    chain_map = {c.id: c for c in chains}
    for pair in all_flat_pairs:
        chain = chain_map[pair.chain_id]
        short_tid = chain.representative_tid[-chain.short_tid_length:]
        if (
            pair.forward_start is not None
            and pair.forward_len is not None
            and pair.reverse_start is not None
            and pair.reverse_len is not None
        ):
            a_start = chain.transcript_offset + pair.forward_start + 1
            a_end = chain.transcript_offset + pair.reverse_start + pair.reverse_len
            pair.pair_name = f"{short_tid}:{a_start}-{a_end}"
        else:
            pair.pair_name = f"{short_tid}:?"

    # ── Step 3: Filter with tnBLAST ────────────────────────────────────────
    log.info("Running tnBLAST specificity filter...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        assay_path = os.path.join(tmp_dir, "assays.txt")
        write_assay_file(all_flat_pairs, assay_path)

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

        # Filter junction-mode pairs (require 0 genomic hits)
        filtered_junction = filter_specific_pairs(
            junction_candidates,
            genome_counts,
            transcriptome_targets,
            target_gene=resolved_gene,
            junction_mode=True,
        )
        # Filter intron-mode pairs (expect 1 genomic hit — the target)
        filtered_intron = filter_specific_pairs(
            intron_candidates,
            genome_counts,
            transcriptome_targets,
            target_gene=resolved_gene,
            junction_mode=False,
        )

        # Merge, deduplicate, prioritise junction over intron
        filtered = dedup_and_prioritize(
            filtered_junction, filtered_intron, max_pairs=num_return,
        )

        # Assign pair numbers
        for i, pair in enumerate(filtered, start=1):
            pair.pair_number = i

        # ── Step 4: Compute genomic fragment lists ─────────────────────────
        for pair in all_flat_pairs:
            name = pair.pair_name

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
                _compute_tnblast_fragments(
                    pair, name, exons, locus,
                    genome_amplicons, transcriptome_amplicons, transcript_exon_map,
                )

    if not filtered:
        log.warning("No primer pairs passed tnBLAST specificity filter.")
    else:
        log.info("%d/%d pairs passed filtering.", len(filtered), len(all_flat_pairs))

    return PipelineResult(
        chains=chains,
        locus=locus,
        filtered_pairs=filtered,
        all_candidates=all_flat_pairs,
        junction_candidates=junction_candidates,
        intron_candidates=intron_candidates,
        warnings=pipeline_warnings,
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
    junction_candidates = result.junction_candidates
    intron_candidates = result.intron_candidates

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
    log.debug("--- Junction-Mode Candidates (includes filtered-out) ---")
    for pair in junction_candidates:
        name = pair.pair_name or "unnamed"
        mark = " [FILTERED IN]" if pair.pair_number is not None else ""
        log.debug(f"  {name}: {pair.forward_seq} / {pair.reverse_seq}  chain={pair.chain_id}{mark}")
        log.debug(f"    P3 F: {pair.primer3_forward_fragments}")
        log.debug(f"    P3 R: {pair.primer3_reverse_fragments}")
        log.debug(f"    tn F: {pair.tnblast_forward_fragments}")
        log.debug(f"    tn R: {pair.tnblast_reverse_fragments}")

    log.debug("")
    log.debug("--- Intron-Mode Candidates (includes filtered-out) ---")
    for pair in intron_candidates:
        name = pair.pair_name or "unnamed"
        mark = " [FILTERED IN]" if pair.pair_number is not None else ""
        log.debug(f"  {name}: {pair.forward_seq} / {pair.reverse_seq}  chain={pair.chain_id}{mark}")
        log.debug(f"    P3 F: {pair.primer3_forward_fragments}")
        log.debug(f"    P3 R: {pair.primer3_reverse_fragments}")
        log.debug(f"    tn F: {pair.tnblast_forward_fragments}")
        log.debug(f"    tn R: {pair.tnblast_reverse_fragments}")

    log.debug("=" * 70)
    log.debug("")
