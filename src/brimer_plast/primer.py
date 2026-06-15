"""Primer3-py integration for Brimer-PLAST.

Wraps primer3-py's ``design_primers`` function with sensible defaults
and returns typed :class:`PrimerPair` objects.

Provides two design modes:

* **Junction mode** (``design_primers``) — at least one primer must span an
  exon-exon junction.  Passes ``SEQUENCE_OVERLAP_JUNCTION_LIST`` to Primer3
  as a soft penalty, then hard-post-filters.
* **Intron mode** (``design_primers_intron_mode``) — forward and reverse
  primers must land on different exons with total intronic separation
  >1000 bp.  No junction constraint.
* **Dual-mode orchestrator** (``design_primers_dual_mode``) — runs both
  modes, returns filtered results.
"""

from __future__ import annotations

from typing import Any

import primer3

from brimer_plast.models import ExonInfo, PrimerPair
from brimer_plast.sequence import exons_in_template_order

# ── Named constants (single source of truth for default values) ──────────────

PRIMER_NUM_RETURN = 10
PRIMER_OPT_SIZE = 20
PRIMER_MIN_SIZE = 18
PRIMER_MAX_SIZE = 25
PRIMER_OPT_TM = 60.0
PRIMER_MIN_TM = 57.0
PRIMER_MAX_TM = 63.0
PRIMER_MIN_GC = 40.0
PRIMER_MAX_GC = 60.0
PRIMER_PRODUCT_MIN = 80
PRIMER_PRODUCT_MAX = 200

# Minimum intronic separation for intron-spanning mode (bp).
DEFAULT_MIN_INTRONIC_SEPARATION = 1000

# How many extra candidates to request from Primer3 per mode so the
# dedup + truncation step has enough to fill the user's --num-return cap.
EXTRA_CANDIDATES = 500


# Default primer3 global args for PCR primer design.
DEFAULT_PRIMER_ARGS: dict[str, Any] = {
    "PRIMER_PRODUCT_SIZE_RANGE": f"{PRIMER_PRODUCT_MIN}-{PRIMER_PRODUCT_MAX}",
    "PRIMER_NUM_RETURN": PRIMER_NUM_RETURN,
    "PRIMER_OPT_SIZE": PRIMER_OPT_SIZE,
    "PRIMER_MIN_SIZE": PRIMER_MIN_SIZE,
    "PRIMER_MAX_SIZE": PRIMER_MAX_SIZE,
    "PRIMER_OPT_TM": PRIMER_OPT_TM,
    "PRIMER_MIN_TM": PRIMER_MIN_TM,
    "PRIMER_MAX_TM": PRIMER_MAX_TM,
    "PRIMER_MIN_GC": PRIMER_MIN_GC,
    "PRIMER_MAX_GC": PRIMER_MAX_GC,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNTP_CONC": 0.6,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MIN_THREE_PRIME_DISTANCE": 3,
}


# ── helpers ──────────────────────────────────────────────────────────────────


def _pair_spans_any_junction(
    raw: dict[str, Any],
    pair_index: int,
    junction_positions_1based: set[int],
) -> bool:
    """Check if at least one primer in the pair spans a required junction.

    Both ``PRIMER_LEFT_*`` and ``PRIMER_RIGHT_*`` coordinates are 0-based
    from the template start.  A primer at 0-based start *s* with length *L*
    spans a 1-based junction *j* iff ``s < j - 1 < s + L``.
    """
    left_start, left_len = raw.get(f"PRIMER_LEFT_{pair_index}", [0, 0])
    right_start, right_len = raw.get(f"PRIMER_RIGHT_{pair_index}", [0, 0])

    for j in junction_positions_1based:
        j_0based = j - 1
        if left_start < j_0based < left_start + left_len:
            return True
        if right_start < j_0based < right_start + right_len:
            return True

    return False


def _build_pair_from_raw(
    raw: dict[str, Any],
    pair_index: int,
    chain_id: str,
) -> PrimerPair | None:
    """Build a :class:`PrimerPair` from primer3-py's raw output.

    Returns ``None`` if the pair data is missing required fields.
    """
    left_pos = raw.get(f"PRIMER_LEFT_{pair_index}")
    right_pos = raw.get(f"PRIMER_RIGHT_{pair_index}")
    if left_pos is None or right_pos is None:
        return None
    return PrimerPair(
        forward_seq=raw.get(f"PRIMER_LEFT_{pair_index}_SEQUENCE", ""),
        reverse_seq=raw.get(f"PRIMER_RIGHT_{pair_index}_SEQUENCE", ""),
        forward_tm=raw.get(f"PRIMER_LEFT_{pair_index}_TM"),
        reverse_tm=raw.get(f"PRIMER_RIGHT_{pair_index}_TM"),
        forward_gc=raw.get(f"PRIMER_LEFT_{pair_index}_GC_PERCENT"),
        reverse_gc=raw.get(f"PRIMER_RIGHT_{pair_index}_GC_PERCENT"),
        product_size=raw.get(f"PRIMER_PAIR_{pair_index}_PRODUCT_SIZE"),
        pair_penalty=raw.get(f"PRIMER_PAIR_{pair_index}_PENALTY"),
        forward_start=left_pos[0],
        forward_len=left_pos[1],
        reverse_start=right_pos[0],
        reverse_len=right_pos[1],
        chain_id=chain_id,
    )


# ── exon-index / intron helpers ──────────────────────────────────────────────


def _exon_index_at_template_pos(
    template_pos_0based: int,
    exons_template_order: list[ExonInfo],
) -> int:
    """Return the template-order index of the exon containing *template_pos_0based*.

    Raises ValueError if the position falls outside the template bounds.
    """
    cumulative = 0
    for i, ex in enumerate(exons_template_order):
        exon_len = ex.end - ex.start + 1
        if cumulative <= template_pos_0based < cumulative + exon_len:
            return i
        cumulative += exon_len
    raise ValueError(
        f"Position {template_pos_0based} is outside the template "
        f"(total length {cumulative} bp)."
    )


def _intronic_separation(
    fwd_start_0based: int,
    rev_start_0based: int,
    exons: list[ExonInfo],
) -> int:
    """Sum of intronic bases between the exons containing forward and reverse primers.

    The forward and reverse primers must reside on *different* exons
    (caller's responsibility to check).  Returns 0 if they end up in the
    same exon, which would fail any >0 threshold.
    """
    ordered = exons_in_template_order(exons)
    # Exons sorted by genomic start (ascending) — works for both strands
    genomic_order = sorted(exons, key=lambda e: e.start)
    n = len(exons)
    strand = exons[0].strand

    def _t_to_g(t_idx: int) -> int:
        """Map template-order index to genomic-order index."""
        if strand == "-":
            return n - 1 - t_idx
        return t_idx

    try:
        fwd_t = _exon_index_at_template_pos(fwd_start_0based, ordered)
        rev_t = _exon_index_at_template_pos(rev_start_0based, ordered)
    except ValueError:
        return 0

    fwd_g = _t_to_g(fwd_t)
    rev_g = _t_to_g(rev_t)

    if fwd_g == rev_g:
        return 0  # same exon

    lo = min(fwd_g, rev_g)
    hi = max(fwd_g, rev_g)

    total = 0
    for i in range(lo, hi):
        cur = genomic_order[i]
        nxt = genomic_order[i + 1]
        total += nxt.start - cur.end - 1

    return total


# ── junction mode (unchanged interface) ──────────────────────────────────────


def design_primers(
    template: str,
    sequence_id: str = "target",
    chain_id: str = "",
    global_args: dict[str, Any] | None = None,
    junction_positions: list[int] | None = None,
    required_junction_positions: list[int] | None = None,
) -> list[PrimerPair]:
    """Design primer pairs for a template sequence using primer3-py.

    Primer3 uses *junction_positions* as a soft penalty via
    ``SEQUENCE_OVERLAP_JUNCTION_LIST``.  This function additionally
    post-filters the results to enforce that at least one primer in each
    pair overlaps at least one required junction (see
    *required_junction_positions*).

    Args:
        template: The spliced exon sequence (DNA string).
        sequence_id: An identifier for the template.
        chain_id: Chain identifier used to populate ``PrimerPair.chain_id``.
        global_args: Overrides for default primer3 global args.
        junction_positions: 1-based positions passed to Primer3 as
            ``SEQUENCE_OVERLAP_JUNCTION_LIST`` (soft penalty).
        required_junction_positions: Subset of *junction_positions* that
            at least one primer must actually overlap for the pair to
            survive post-filtering.  When ``None``, defaults to
            *junction_positions* (i.e. all listed junctions are required).

    Returns:
        A list of :class:`PrimerPair` objects, one per candidate pair.
        Empty list if no primers could be designed or none passed
        the junction-overlap post-filter.
    """
    merged_args = {**DEFAULT_PRIMER_ARGS, **(global_args or {})}

    seq_args: dict[str, Any] = {
        "SEQUENCE_ID": sequence_id,
        "SEQUENCE_TEMPLATE": template,
    }

    if junction_positions:
        seq_args["SEQUENCE_OVERLAP_JUNCTION_LIST"] = junction_positions
        merged_args.setdefault("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", 1)

    required_set: set[int] | None = None
    if required_junction_positions is not None:
        required_set = set(required_junction_positions)
    elif junction_positions:
        required_set = set(junction_positions)

    try:
        raw: dict[str, Any] = primer3.design_primers(seq_args, merged_args)
    except OSError:
        return []

    num_pairs = raw.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if not num_pairs:
        return []

    result: list[PrimerPair] = []
    for i in range(num_pairs):
        pair = _build_pair_from_raw(raw, i, chain_id)
        if pair is None:
            continue
        if required_set is not None and not _pair_spans_any_junction(raw, i, required_set):
            continue
        result.append(pair)

    return result


# ── intron mode ──────────────────────────────────────────────────────────────


def design_primers_intron_mode(
    template: str,
    exons: list[ExonInfo],
    sequence_id: str = "target",
    chain_id: str = "",
    global_args: dict[str, Any] | None = None,
    min_intronic_separation: int = DEFAULT_MIN_INTRONIC_SEPARATION,
) -> list[PrimerPair]:
    """Design primer pairs separated by introns > *min_intronic_separation*.

    No ``SEQUENCE_OVERLAP_JUNCTION_LIST`` is passed to Primer3 — the only
    constraint is that forward and reverse primers fall on different exons
    and the sum of intronic bases between them exceeds the threshold.

    Args:
        template: The spliced exon sequence (DNA string).
        exons: Exon list for this chain (used for coordinate mapping).
        sequence_id: An identifier for the template.
        chain_id: Chain identifier used to populate ``PrimerPair.chain_id``.
        global_args: Overrides for default primer3 global args.
        min_intronic_separation: Minimum total intronic bases (default 1000).

    Returns:
        A list of :class:`PrimerPair` objects, one per passing pair.
        Empty list if no primers could be designed or none passed.
    """
    if len(exons) < 2:
        return []  # single-exon chain can never have intronic separation

    merged_args = {**DEFAULT_PRIMER_ARGS, **(global_args or {})}

    seq_args: dict[str, Any] = {
        "SEQUENCE_ID": sequence_id,
        "SEQUENCE_TEMPLATE": template,
    }

    try:
        raw: dict[str, Any] = primer3.design_primers(seq_args, merged_args)
    except OSError:
        return []

    num_pairs = raw.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if not num_pairs:
        return []

    result: list[PrimerPair] = []
    for i in range(num_pairs):
        pair = _build_pair_from_raw(raw, i, chain_id)
        if pair is None:
            continue
        fwd_start = pair.forward_start
        rev_start = pair.reverse_start
        if fwd_start is None or rev_start is None:
            continue

        # Must be on different exons
        ordered = exons_in_template_order(exons)
        try:
            fwd_idx = _exon_index_at_template_pos(fwd_start, ordered)
            rev_idx = _exon_index_at_template_pos(rev_start, ordered)
        except ValueError:
            continue

        if fwd_idx == rev_idx:
            continue

        separation = _intronic_separation(fwd_start, rev_start, exons)
        if separation <= min_intronic_separation:
            continue

        result.append(pair)

    return result


# ── dual-mode orchestrator ───────────────────────────────────────────────────


def design_primers_dual_mode(
    template: str,
    exons: list[ExonInfo],
    sequence_id: str = "target",
    chain_id: str = "",
    global_args: dict[str, Any] | None = None,
    num_return: int = PRIMER_NUM_RETURN,
    junction_positions: list[int] | None = None,
    required_junction_positions: list[int] | None = None,
    min_intronic_separation: int = DEFAULT_MIN_INTRONIC_SEPARATION,
) -> tuple[list[PrimerPair], list[PrimerPair]]:
    """Design primers in both junction and intron modes.

    Both modes request ``num_return + EXTRA_CANDIDATES`` from Primer3 to
    ensure the downstream dedup step has enough candidates to fill the
    user's requested cap.

    Returns:
        ``(junction_pairs, intron_pairs)`` — pre-filtered lists (each
        pair has passed its respective mode's post-filter).
    """
    boosted_args = {**(global_args or {})}
    boosted_args["PRIMER_NUM_RETURN"] = num_return + EXTRA_CANDIDATES

    # Mode A — junction-spanning
    junction_pairs = design_primers(
        template,
        sequence_id=sequence_id,
        chain_id=chain_id,
        global_args=boosted_args,
        junction_positions=junction_positions,
        required_junction_positions=required_junction_positions,
    )

    # Mode B — intron-spanning
    intron_pairs = design_primers_intron_mode(
        template,
        exons,
        sequence_id=sequence_id,
        chain_id=chain_id,
        global_args=boosted_args,
        min_intronic_separation=min_intronic_separation,
    )

    return junction_pairs, intron_pairs


# ── post-pipeline dedup ──────────────────────────────────────────────────────


def dedup_and_prioritize(
    junction_pairs: list[PrimerPair],
    intron_pairs: list[PrimerPair],
    max_pairs: int = PRIMER_NUM_RETURN,
) -> list[PrimerPair]:
    """Merge two lists, deduplicate by primer sequence, prioritise junction.

    Junction-mode pairs come first in the merged list, so a pair that
    passes both filters appears under its junction-mode entry.  The
    dedup key is ``(forward_seq, reverse_seq)``.

    Args:
        junction_pairs: Pairs from junction mode.
        intron_pairs: Pairs from intron mode.
        max_pairs: Maximum pairs to return (output cap).

    Returns:
        Deduped, prioritised list, truncated at *max_pairs*.
    """
    seen: set[tuple[str, str]] = set()
    result: list[PrimerPair] = []

    for pair in junction_pairs + intron_pairs:
        key = (pair.forward_seq, pair.reverse_seq)
        if key in seen:
            continue
        seen.add(key)
        result.append(pair)
        if len(result) >= max_pairs:
            break

    return result
