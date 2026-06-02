"""Primer3-py integration for Brimer-PLAST.

Wraps primer3-py's ``design_primers`` function with sensible defaults
and returns typed :class:`PrimerPair` objects.
"""

from typing import Any

import primer3

from brimer_plast.models import PrimerPair


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


# Default primer3 global args for PCR primer design.
DEFAULT_PRIMER_ARGS: dict[str, Any] = {
    "PRIMER_PRODUCT_SIZE_RANGE": "100-400",
    "PRIMER_NUM_RETURN": 50,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNTP_CONC": 0.6,
    "PRIMER_DNA_CONC": 50.0,
}


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
        left_pos = raw.get(f"PRIMER_LEFT_{i}", [0, 0])
        right_pos = raw.get(f"PRIMER_RIGHT_{i}", [0, 0])
        pair = PrimerPair(
            forward_seq=raw.get(f"PRIMER_LEFT_{i}_SEQUENCE", ""),
            reverse_seq=raw.get(f"PRIMER_RIGHT_{i}_SEQUENCE", ""),
            forward_tm=raw.get(f"PRIMER_LEFT_{i}_TM"),
            reverse_tm=raw.get(f"PRIMER_RIGHT_{i}_TM"),
            forward_gc=raw.get(f"PRIMER_LEFT_{i}_GC_PERCENT"),
            reverse_gc=raw.get(f"PRIMER_RIGHT_{i}_GC_PERCENT"),
            product_size=raw.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
            pair_penalty=raw.get(f"PRIMER_PAIR_{i}_PENALTY"),
            forward_start=left_pos[0],
            forward_len=left_pos[1],
            reverse_start=right_pos[0],
            reverse_len=right_pos[1],
            chain_id=chain_id,
        )
        if required_set is not None and not _pair_spans_any_junction(
            raw, i, required_set
        ):
            continue
        result.append(pair)

    return result
