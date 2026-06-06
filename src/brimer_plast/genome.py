"""Conserved exon chain detection and top-level orchestration for Brimer-PLAST.

This module is the main entry point for genome-related operations.  It
provides:

* ``get_target_information`` — the primary entry point for extracting
  conserved exon chains from a genome + annotation for a given target.
* ``compute_conserved_exon_chains`` — conserved chain detection logic.
* ``get_gene_locus`` — builds the full gene locus for report visualisation.

Lower-level GTF parsing and sequence extraction live in ``gtf.py`` and
``sequence.py`` respectively, but are re-exported here for backward
compatibility.
"""

from __future__ import annotations

from pathlib import Path

import pyfaidx

from brimer_plast.gtf import (
    _find_gene_for_transcript,
    build_transcript_to_gene_map,
    parse_gtf,
    parse_gtf_all_transcripts,
    parse_gtf_grouped_by_transcript,
)
from brimer_plast.models import ConservedExonChain, ExonInfo, GeneLocus, GenomicFragment
from brimer_plast.sequence import (
    _extract_sequence_from_genome,
    build_transcriptome_fasta,
    exons_in_template_order,
    extract_sequence,
    genomic_range_to_fragments,
    reverse_complement,
    template_to_genomic,
)

# ── Re-exports for backward compatibility ────────────────────────────────────
__all__ = [
    # from gtf
    "build_transcript_to_gene_map",
    "parse_gtf",
    "parse_gtf_all_transcripts",
    "parse_gtf_grouped_by_transcript",
    # from sequence
    "build_transcriptome_fasta",
    "exons_in_template_order",
    "extract_sequence",
    "genomic_range_to_fragments",
    "reverse_complement",
    "template_to_genomic",
    # own
    "compute_conserved_exon_chains",
    "get_gene_locus",
    "get_target_information",
]


# ── exon key ─────────────────────────────────────────────────────────────────


def _exon_key(exon: ExonInfo) -> tuple[int, int]:
    """Return a stable identifier for an exon: (start, end)."""
    return (exon.start, exon.end)


# ── junction computation helpers ─────────────────────────────────────────────


def _compute_junction_positions(exons: list[ExonInfo]) -> list[int]:
    """Return 1-based junction positions within the template.

    Assumes *exons* are already in template order.
    """
    if len(exons) < 2:
        return []
    positions: list[int] = []
    cumulative_len = 0
    for i in range(len(exons) - 1):
        ex = exons[i]
        exon_len = ex.end - ex.start + 1
        cumulative_len += exon_len
        positions.append(cumulative_len + 1)
    return positions


def _compute_junction_adjacencies(
    exons: list[ExonInfo],
) -> set[tuple[tuple[int, int], tuple[int, int]]]:
    """Return the set of exon-key adjacency pairs for a list of exons.

    Exons should be in template order.  Each adjacency is
    ``((exon_i_start, exon_i_end), (exon_{i+1}_start, exon_{i+1}_end))``.
    """
    result: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for i in range(len(exons) - 1):
        result.add((_exon_key(exons[i]), _exon_key(exons[i + 1])))
    return result


def _compute_unique_junction_positions(
    all_junctions: list[int],
    target_exons: list[ExonInfo],
    sibling_exon_lists: list[list[ExonInfo]],
) -> list[int]:
    """Return 1-based junction positions unique to the target transcript.

    Only junctions (adjacent exon pairs that appear in the target transcript
    but NOT in any sibling transcript) are kept.  *target_exons* must be in
    template order.  If no unique junctions exist, returns an empty list.
    """
    target_adj = _compute_junction_adjacencies(target_exons)

    sibling_adj: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for sib in sibling_exon_lists:
        sibling_adj |= _compute_junction_adjacencies(exons_in_template_order(sib))

    unique_adj = target_adj - sibling_adj
    if not unique_adj:
        return []

    # Map unique adjacencies back to 1-based template positions
    # by walking the target exons in template order
    result: list[int] = []
    cumulative_len = 0
    for i in range(len(target_exons) - 1):
        exon_len = target_exons[i].end - target_exons[i].start + 1
        cumulative_len += exon_len
        adj = (_exon_key(target_exons[i]), _exon_key(target_exons[i + 1]))
        if adj in unique_adj:
            result.append(cumulative_len + 1)

    return result


# ── conserved exon chains ────────────────────────────────────────────────────


def compute_conserved_exon_chains(
    transcript_exon_lists: list[list[ExonInfo]],
) -> list[ConservedExonChain]:
    """Find contiguous exon chains conserved across all transcripts.

    A *contiguous conserved chain* is a maximal run of exons (in genomic
    order) where every adjacent pair's *exon-key adjacency* appears in all
    transcripts.  Exon-key = ``(start, end)``.

    Args:
        transcript_exon_lists: exon lists for each transcript (must share
            the same seqid and strand).

    Returns:
        List of :class:`ConservedExonChain` objects.  Template and
        junction positions are NOT yet populated (call
        :func:`get_target_information` for that).

    Raises:
        ValueError: if input is empty, transcripts span different seqids
            or strands, any transcript has < 2 exons, or no conserved
            adjacencies exist.
    """
    if not transcript_exon_lists:
        raise ValueError("No transcript exon lists provided.")

    # Validate shared seqid / strand
    first = transcript_exon_lists[0]
    if not first:
        raise ValueError("Empty exon list for first transcript.")
    ref_seqid = first[0].seqid
    ref_strand = first[0].strand
    for tlist in transcript_exon_lists:
        for ex in tlist:
            if ex.seqid != ref_seqid:
                raise ValueError(f"Transcripts span different chromosomes: {ex.seqid!r}")
            if ex.strand != ref_strand:
                raise ValueError(f"Transcripts on different strands: {ex.strand!r}")

    # Reject single-exon transcripts (no possible junctions)
    for tidx, tlist in enumerate(transcript_exon_lists):
        if len(tlist) < 2:
            raise ValueError(
                f"Transcript at index {tidx} has {len(tlist)} exon(s); "
                "at least 2 exons are needed for an exon-exon junction. "
                "Use --disable-junction-overlap to design primers "
                "without junction-spanning."
            )

    # Build adjacency sets using (start, end) as exon key
    all_adjacencies: list[set[tuple[tuple[int, int], tuple[int, int]]]] = []
    for tlist in transcript_exon_lists:
        adj: set[tuple[tuple[int, int], tuple[int, int]]] = set()
        for i in range(len(tlist) - 1):
            adj.add((_exon_key(tlist[i]), _exon_key(tlist[i + 1])))
        all_adjacencies.append(adj)

    conserved_adj = all_adjacencies[0]
    for adj in all_adjacencies[1:]:
        conserved_adj &= adj

    if not conserved_adj:
        raise ValueError(
            "No conserved exon-exon junctions found across transcripts. "
            "Use --disable-junction-overlap to design primers "
            "without junction-spanning."
        )

    # Build a lookup from exon key -> ExonInfo (union across transcripts)
    exon_map: dict[tuple[int, int], ExonInfo] = {}
    for tlist in transcript_exon_lists:
        for ex in tlist:
            exon_map[_exon_key(ex)] = ex

    # Build successor map: left_key -> right_key
    succ: dict[tuple[int, int], tuple[int, int]] = {}
    for lk, rk in conserved_adj:
        succ[lk] = rk

    # Chain starts: left keys that are not also right keys
    left_keys = {a[0] for a in conserved_adj}
    right_keys = {a[1] for a in conserved_adj}
    chain_starts = left_keys - right_keys

    if not chain_starts:
        raise ValueError(
            "No conserved exon-exon junctions found across transcripts. "
            "Use --disable-junction-overlap to design primers "
            "without junction-spanning."
        )

    chains: list[ConservedExonChain] = []

    for start_key in sorted(chain_starts):
        exons: list[ExonInfo] = []
        cur_key: tuple[int, int] | None = start_key
        while cur_key is not None:
            ex = exon_map.get(cur_key)
            if ex is None:
                break
            exons.append(ex)
            cur_key = succ.get(cur_key)

        if len(exons) >= 2:
            chains.append(
                ConservedExonChain(
                    id="",
                    exons=exons,
                    template="",
                    junction_positions_1based=[],
                )
            )

    if not chains:
        raise ValueError(
            "No conserved exon-exon junctions found across transcripts. "
            "Use --disable-junction-overlap to design primers "
            "without junction-spanning."
        )

    return chains


# ── gene locus ───────────────────────────────────────────────────────────────


def get_gene_locus(
    gtf_path: str | Path,
    target_gene: str | None = None,
    target_transcript: str | None = None,
) -> GeneLocus | None:
    """Fetch all transcripts and coordinates for the target gene/transcript.

    If target_transcript is provided, it finds the sibling gene and returns
    all transcripts for that gene.
    """
    resolved_gene = target_gene
    if target_transcript:
        resolved_gene = _find_gene_for_transcript(gtf_path, target_transcript)

    if not resolved_gene:
        return None

    transcripts = parse_gtf_grouped_by_transcript(gtf_path, target_gene=resolved_gene)
    if not transcripts:
        return None

    # Deduplicate transcripts by exon structure (same seqid, strand, and starts/ends)
    unique_transcripts: dict[str, list[ExonInfo]] = {}
    seen_structures = set()
    for tid, exons in transcripts.items():
        struct = tuple((e.start, e.end) for e in exons)
        if struct not in seen_structures:
            unique_transcripts[tid] = exons
            seen_structures.add(struct)

    first_exons = list(unique_transcripts.values())[0]
    seqid = first_exons[0].seqid
    strand = first_exons[0].strand

    all_starts = []
    all_ends = []
    for exons in unique_transcripts.values():
        for e in exons:
            all_starts.append(e.start)
            all_ends.append(e.end)

    return GeneLocus(
        gene_name=resolved_gene,
        seqid=seqid,
        strand=strand,
        transcripts=unique_transcripts,
        min_start=min(all_starts),
        max_end=max(all_ends),
    )


# ── helper: transcript_id lookup ────────────────────────────────────────────


def _tid_for_exons(
    transcripts: dict[str, list[ExonInfo]],
    exons: list[ExonInfo],
) -> str:
    """Look up the transcript_id for an exon list."""
    for tid, el in transcripts.items():
        if el is exons or el == exons:
            return tid
    return ""


# ── naming helpers ─────────────────────────────────────────────────────────


def _compute_short_tid_length(transcript_ids: list[str]) -> int:
    """Fifth-shortest suffix length that uniquely identifies every ID.

    Finds the smallest L such that the last L characters of every
    transcript ID are unique within the list, then returns L + 4.
    With a single transcript, L=1 trivially, so returns 5.
    """
    if not transcript_ids:
        return 0
    if len(transcript_ids) == 1:
        return 5

    max_len = max(len(t) for t in transcript_ids)
    for L in range(1, max_len + 1):
        suffixes = [t[-L:] for t in transcript_ids]
        if len(set(suffixes)) == len(transcript_ids):
            return L + 4
    return max_len + 4


def _transcript_contains_chain(
    transcript_exons: list[ExonInfo],
    chain_exons: list[ExonInfo],
) -> bool:
    """True if every exon in *chain_exons* (by start/end key) is present."""
    transcript_keys = {_exon_key(ex) for ex in transcript_exons}
    return all(_exon_key(ex) in transcript_keys for ex in chain_exons)


def _compute_chain_offset_in_transcript(
    transcript_exons: list[ExonInfo],
    chain_first_exon: ExonInfo,
) -> int:
    """Cumulative bp of transcript exons (5'\u21923') before chain's first exon."""
    ordered = exons_in_template_order(transcript_exons)
    offset = 0
    for ex in ordered:
        if _exon_key(ex) == _exon_key(chain_first_exon):
            break
        offset += ex.end - ex.start + 1
    else:
        raise ValueError(
            f"Chain's first exon ({_exon_key(chain_first_exon)}) "
            f"not found in the chosen transcript."
        )
    return offset


def _pick_representative_for_chain(
    chain: ConservedExonChain,
    transcripts: dict[str, list[ExonInfo]],
) -> tuple[str, int]:
    """Pick the alphanumerically first transcript that contains *chain*, and
    return ``(representative_tid, transcript_offset)``.

    Raises ValueError if no transcript contains the chain's exons.
    """
    candidates = [
        tid
        for tid, exons in transcripts.items()
        if _transcript_contains_chain(exons, chain.exons)
    ]
    if not candidates:
        raise ValueError(
            f"No transcript contains chain {chain.id!r}. This should "
            f"not happen — the chain was derived from these transcripts."
        )

    rep = min(candidates)  # alphanumerically first
    offset = _compute_chain_offset_in_transcript(
        transcripts[rep], chain.exons[0],
    )
    return rep, offset


# ── per-mode helpers ─────────────────────────────────────────────────────────


def _build_target_transcript_chain(
    genome: pyfaidx.Fasta,
    target_transcript: str,
    gtf_path: str | Path,
    transcripts: dict[str, list[ExonInfo]],
    transcript_exon_lists: list[list[ExonInfo]],
) -> list[ConservedExonChain]:
    """Build a single conserved exon chain for a target transcript.

    All junctions are passed to Primer3, but only junctions *unique* to
    this transcript (not shared with any sibling of the same gene) are
    marked as required.
    """
    exons = transcript_exon_lists[0]
    template = _extract_sequence_from_genome(genome, exons)

    if len(exons) < 2:
        return [
            ConservedExonChain(
                id=target_transcript,
                exons=exons,
                template=template,
                junction_positions_1based=[],
                required_junction_positions_1based=[],
            )
        ]

    template_order_exons = exons_in_template_order(exons)
    junctions = _compute_junction_positions(template_order_exons)

    gene_name = _find_gene_for_transcript(gtf_path, target_transcript)
    if gene_name:
        all_siblings = parse_gtf_grouped_by_transcript(gtf_path, target_gene=gene_name)
        sibling_exon_lists = [
            t for tid, t in all_siblings.items() if tid != target_transcript
        ]
        if sibling_exon_lists:
            unique_positions = _compute_unique_junction_positions(
                junctions, template_order_exons, sibling_exon_lists
            )
            required_junctions = unique_positions
        else:
            required_junctions = junctions
    else:
        required_junctions = junctions

    return [
        ConservedExonChain(
            id=target_transcript,
            exons=exons,
            template=template,
            junction_positions_1based=junctions,
            required_junction_positions_1based=required_junctions,
            transcript_offset=0,
            representative_tid=target_transcript,
            short_tid_length=5,
        )
    ]


def _build_target_gene_chains(
    genome: pyfaidx.Fasta,
    target_gene: str,
    transcripts: dict[str, list[ExonInfo]],
    transcript_exon_lists: list[list[ExonInfo]],
) -> list[ConservedExonChain]:
    """Build conserved exon chains for a target gene.

    Multi-exon transcripts are analysed for conserved exon chains.
    If conserved chains exist, those are returned along with any
    single-exon transcripts as junctionless chains.  If no conserved
    adjacencies are found, each multi-exon transcript becomes its own chain.

    Each chain is annotated with transcript_offset, representative_tid,
    and short_tid_length for primer pair naming (see ADR 0002).
    """
    multi_exon_lists = [tl for tl in transcript_exon_lists if len(tl) >= 2]
    single_exon_lists = [tl for tl in transcript_exon_lists if len(tl) < 2]

    all_tids = list(transcripts.keys())
    short_tid_length = _compute_short_tid_length(all_tids)

    result: list[ConservedExonChain] = []

    if multi_exon_lists:
        try:
            chains = compute_conserved_exon_chains(multi_exon_lists)
            for chain_idx, chain in enumerate(chains, start=1):
                template = _extract_sequence_from_genome(genome, chain.exons)
                template_order_exons = exons_in_template_order(chain.exons)
                junctions = _compute_junction_positions(template_order_exons)
                chain_id = f"{target_gene}_chain_{chain_idx}"
                rep, offset = _pick_representative_for_chain(
                    chain, transcripts
                )
                result.append(
                    ConservedExonChain(
                        id=chain_id,
                        exons=chain.exons,
                        template=template,
                        junction_positions_1based=junctions,
                        required_junction_positions_1based=junctions,
                        transcript_offset=offset,
                        representative_tid=rep,
                        short_tid_length=short_tid_length,
                    )
                )
        except ValueError:
            for tl in multi_exon_lists:
                template = _extract_sequence_from_genome(genome, tl)
                template_order_exons = exons_in_template_order(tl)
                junctions = _compute_junction_positions(template_order_exons)
                tid = _tid_for_exons(transcripts, tl)
                result.append(
                    ConservedExonChain(
                        id=tid,
                        exons=tl,
                        template=template,
                        junction_positions_1based=junctions,
                        required_junction_positions_1based=junctions,
                        transcript_offset=0,
                        representative_tid=tid,
                        short_tid_length=short_tid_length,
                        fallback=True,
                    )
                )

    for tl in single_exon_lists:
        template = _extract_sequence_from_genome(genome, tl)
        tid = _tid_for_exons(transcripts, tl)
        result.append(
            ConservedExonChain(
                id=tid,
                exons=tl,
                template=template,
                junction_positions_1based=[],
                required_junction_positions_1based=[],
                transcript_offset=0,
                representative_tid=tid,
                short_tid_length=short_tid_length,
            )
        )

    return result


# ── top-level entry point ────────────────────────────────────────────────────


def get_target_information(
    fasta_path: str | Path,
    gtf_path: str | Path,
    target_gene: str | None = None,
    target_transcript: str | None = None,
) -> list[ConservedExonChain]:
    """Get conserved exon chains with spliced template sequences.

    Exactly one of *target_gene* or *target_transcript* must be provided.

    Delegates to :func:`_build_target_transcript_chain` or
    :func:`_build_target_gene_chains` based on the target type.

    Returns a list of :class:`ConservedExonChain` objects with populated
    templates and junction positions.
    """
    if target_gene and target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript, not both.")
    if not target_gene and not target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript.")

    transcripts = parse_gtf_grouped_by_transcript(
        gtf_path,
        target_gene=target_gene,
        target_transcript=target_transcript,
    )
    transcript_exon_lists = list(transcripts.values())

    genome = pyfaidx.Fasta(str(fasta_path), read_ahead=10_000)
    try:
        if target_transcript:
            return _build_target_transcript_chain(
                genome, target_transcript, gtf_path, transcripts, transcript_exon_lists,
            )
        return _build_target_gene_chains(
            genome, target_gene or "", transcripts, transcript_exon_lists,
        )
    finally:
        genome.close()