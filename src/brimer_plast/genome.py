"""GTF parsing and genome sequence extraction for Brimer-PLAST."""

import csv
import re
from pathlib import Path

import pyfaidx

from brimer_plast.models import ConservedExonChain, ExonInfo


# ── flat parsing (backward-compat) ───────────────────────────────────────────


def parse_gtf(
    gtf_path: str | Path,
    target_gene: str | None = None,
    target_transcript: str | None = None,
) -> list[ExonInfo]:
    """Parse a GTF file and return exon coordinates for the target.

    Exactly one of *target_gene* or *target_transcript* must be provided.

    Returns a flat list of :class:`ExonInfo` objects with 1-indexed coordinates,
    matching all exons across all transcripts (when given a gene name).
    """
    if target_gene and target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript, not both.")
    if not target_gene and not target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript.")

    exons: list[ExonInfo] = []

    with open(gtf_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 9:
                continue
            seqid, source, feature_type, start_s, end_s, _, strand, _, attr = row

            if feature_type != "exon":
                continue

            start = int(start_s)
            end = int(end_s)
            attrs = _parse_gtf_attributes(attr)

            if target_gene:
                gene_val = attrs.get("gene_name") or attrs.get("gene") or attrs.get("gene_id")
                if gene_val == target_gene:
                    exons.append(ExonInfo(seqid=seqid, start=start, end=end, strand=strand))
            elif target_transcript:
                if attrs.get("transcript_id") == target_transcript:
                    exons.append(ExonInfo(seqid=seqid, start=start, end=end, strand=strand))

    if not exons:
        key = target_gene or target_transcript or ""
        raise ValueError(f"Target {key!r} not found in {gtf_path}")

    return exons


# ── transcript-grouped parsing ───────────────────────────────────────────────


def parse_gtf_grouped_by_transcript(
    gtf_path: str | Path,
    target_gene: str | None = None,
    target_transcript: str | None = None,
) -> dict[str, list[ExonInfo]]:
    """Parse a GTF and return exons grouped by ``transcript_id``.

    Exactly one of *target_gene* or *target_transcript* must be provided.

    Returns ``{transcript_id: [ExonInfo, ...], ...}``.  When given a gene
    name, ALL transcripts of that gene are returned.  When given a transcript
    ID, the dict contains a single entry.
    """
    if target_gene and target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript, not both.")
    if not target_gene and not target_transcript:
        raise ValueError("Provide either --target-gene or --target-transcript.")

    result: dict[str, list[ExonInfo]] = {}

    with open(gtf_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 9:
                continue
            seqid, source, feature_type, start_s, end_s, _, strand, _, attr = row

            if feature_type != "exon":
                continue

            start = int(start_s)
            end = int(end_s)
            attrs = _parse_gtf_attributes(attr)
            tid = attrs.get("transcript_id")
            if not tid:
                continue

            if target_gene:
                gene_val = attrs.get("gene_name") or attrs.get("gene") or attrs.get("gene_id")
                if gene_val != target_gene:
                    continue
            elif target_transcript:
                if tid != target_transcript:
                    continue

            result.setdefault(tid, []).append(
                ExonInfo(seqid=seqid, start=start, end=end, strand=strand)
            )

    if not result:
        key = target_gene or target_transcript or ""
        raise ValueError(f"Target {key!r} not found in {gtf_path}")

    for tid in result:
        result[tid].sort(key=lambda e: e.start)

    return result


# ── conserved exon chains ────────────────────────────────────────────────────


def _exon_key(exon: ExonInfo) -> tuple[int, int]:
    """Return a stable identifier for an exon: (start, end)."""
    return (exon.start, exon.end)


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


# ── template extraction ──────────────────────────────────────────────────────


def extract_sequence(
    fasta_path: str | Path,
    exons: list[ExonInfo],
) -> str:
    """Extract and concatenate exon sequences from a FASTA genome.

    For exons on the negative strand, the sequence is reverse-complemented.
    Exon coordinates are 1-indexed (GTF convention).

    Raises:
        ValueError: If exons span different chromosomes, or seqid not found.
    """
    if not exons:
        return ""

    genome = pyfaidx.Fasta(str(fasta_path), read_ahead=10_000)

    seqid = exons[0].seqid
    for exon in exons:
        if exon.seqid != seqid:
            genome.close()
            raise ValueError(
                f"Exons span different chromosomes: {exon.seqid!r} != {seqid!r}. "
                "All exons for a target must share the same seqid."
            )

    if seqid not in genome:
        raise ValueError(
            f"Sequence {seqid!r} not found in {fasta_path}. "
            f"Available sequences: {list(genome.keys())[:10]}"
        )

    ordered = _exons_in_template_order(exons)

    parts: list[str] = []
    for exon in ordered:
        fragment = genome[seqid][exon.start - 1 : exon.end].seq
        if exon.strand == "-":
            fragment = _reverse_complement(fragment)
        parts.append(fragment)

    genome.close()
    return "".join(parts)


def parse_gtf_all_transcripts(
    gtf_path: str | Path,
) -> dict[str, list[ExonInfo]]:
    """Parse a GTF and return ALL transcripts with their exon lists.

    Returns ``{transcript_id: [ExonInfo, ...], ...}`` for every transcript
    found in the GTF.  Use this for building the full transcriptome.
    """
    result: dict[str, list[ExonInfo]] = {}

    with open(gtf_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 9:
                continue
            seqid, source, feature_type, start_s, end_s, _, strand, _, attr = row

            if feature_type != "exon":
                continue

            start = int(start_s)
            end = int(end_s)
            attrs = _parse_gtf_attributes(attr)
            tid = attrs.get("transcript_id")
            if not tid:
                continue

            result.setdefault(tid, []).append(
                ExonInfo(seqid=seqid, start=start, end=end, strand=strand)
            )

    for tid in result:
        result[tid].sort(key=lambda e: e.start)

    return result


def build_transcriptome_fasta(
    fasta_path: str | Path,
    gtf_path: str | Path,
    output_path: str | Path,
) -> None:
    """Build a transcriptome FASTA by splicing all transcripts from GTF.

    Extracts exon sequences from the genome FASTA and splices them per
    transcript (reverse-complementing exons on the negative strand).
    Writes a FASTA file with one entry per transcript_id.
    """
    transcripts = parse_gtf_all_transcripts(gtf_path)
    with open(output_path, "w") as f:
        for tid, exons in transcripts.items():
            template = extract_sequence(fasta_path, exons)
            f.write(f">{tid}\n")
            # Write 60-char lines for readability
            for i in range(0, len(template), 60):
                f.write(template[i : i + 60] + "\n")


def build_transcript_to_gene_map(
    gtf_path: str | Path,
) -> dict[str, str]:
    """Build a map from transcript_id to gene_name from a GTF.

    Returns ``{transcript_id: gene_name, ...}``.  Falls back to ``gene_id``
    if ``gene_name`` attribute is absent.
    """
    result: dict[str, str] = {}

    with open(gtf_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 9:
                continue
            feature_type = row[2]
            if feature_type != "exon":
                continue
            attrs = _parse_gtf_attributes(row[8])
            tid = attrs.get("transcript_id")
            if not tid or tid in result:
                continue
            gene = attrs.get("gene_name") or attrs.get("gene") or attrs.get("gene_id")
            if gene:
                result[tid] = gene

    return result


def _exons_in_template_order(exons: list[ExonInfo]) -> list[ExonInfo]:
    """Return exons sorted in the order they appear in the spliced template."""
    if not exons:
        return []
    strand = exons[0].strand
    if strand == "-":
        return sorted(exons, key=lambda e: e.start, reverse=True)
    else:
        return sorted(exons, key=lambda e: e.start)


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


# ── top-level entry point ────────────────────────────────────────────────────


def _find_gene_for_transcript(
    gtf_path: str | Path,
    transcript_id: str,
) -> str | None:
    """Scan GTF to find which gene a transcript belongs to.

    Returns the gene name (``gene_name``, ``gene``, or ``gene_id`` from
    GTF attributes) or ``None`` if the transcript is not found or has no
    gene information.
    """
    with open(gtf_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 9:
                continue
            if row[2] != "exon":
                continue
            attrs = _parse_gtf_attributes(row[8])
            if attrs.get("transcript_id") == transcript_id:
                return (attrs.get("gene_name") or attrs.get("gene") or attrs.get("gene_id"))
    return None


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
        sibling_adj |= _compute_junction_adjacencies(
            _exons_in_template_order(sib)
        )

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


def get_target_information(
    fasta_path: str | Path,
    gtf_path: str | Path,
    target_gene: str | None = None,
    target_transcript: str | None = None,
) -> list[ConservedExonChain]:
    """Get conserved exon chains with spliced template sequences.

    Exactly one of *target_gene* or *target_transcript* must be provided.

    This function always returns whatever biological data is available.
    Single-exon transcripts produce junctionless chains.  Transcripts with
    no unique junctions produce chains with empty ``required_junction_positions_1based``.
    The caller (typically the CLI) is responsible for enforcing junction-spanning
    policy (e.g. dropping junctionless chains when --disable-junction-overlap is
    not set).

    When *target_gene* is provided, multi-exon transcripts are analysed for
    conserved exon chains.  If conserved chains exist, those are returned along
    with any single-exon transcripts as junctionless chains.  If no conserved
    adjacencies are found, each multi-exon transcript becomes its own chain.

    When *target_transcript* is provided, returns a single chain for that
    transcript.  All junctions are passed to Primer3, but only junctions
    *unique* to this transcript (not shared with any sibling of the same
    gene) are marked as required.

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

    def _tid_for_exons(exons: list[ExonInfo]) -> str:
        """Look up the transcript_id for an exon list."""
        for tid, el in transcripts.items():
            if el is exons or el == exons:
                return tid
        return ""

    if target_transcript:
        exons = transcript_exon_lists[0]
        template = extract_sequence(fasta_path, exons)

        if len(exons) < 2:
            # Single-exon transcript — junctionless chain.
            return [
                ConservedExonChain(
                    id=target_transcript,
                    exons=exons,
                    template=template,
                    junction_positions_1based=[],
                    required_junction_positions_1based=[],
                )
            ]

        template_order_exons = _exons_in_template_order(exons)
        junctions = _compute_junction_positions(template_order_exons)

        # Try to find unique junctions for transcript-specificity
        gene_name = _find_gene_for_transcript(gtf_path, target_transcript)
        if gene_name:
            all_siblings = parse_gtf_grouped_by_transcript(
                gtf_path, target_gene=gene_name
            )
            sibling_exon_lists = [
                t for tid, t in all_siblings.items() if tid != target_transcript
            ]
            if sibling_exon_lists:
                unique_positions = _compute_unique_junction_positions(
                    junctions, template_order_exons, sibling_exon_lists
                )
                required_junctions = unique_positions  # may be empty
            else:
                # Only transcript of this gene — all junctions are required
                required_junctions = junctions
        else:
            # Fall back: cannot determine gene for this transcript
            required_junctions = junctions

        return [
            ConservedExonChain(
                id=target_transcript,
                exons=exons,
                template=template,
                junction_positions_1based=junctions,
                required_junction_positions_1based=required_junctions,
            )
        ]
    else:
        # Separate single-exon from multi-exon transcripts
        multi_exon_lists = [
            tl for tl in transcript_exon_lists if len(tl) >= 2
        ]
        single_exon_lists = [
            tl for tl in transcript_exon_lists if len(tl) < 2
        ]

        result: list[ConservedExonChain] = []

        if multi_exon_lists:
            try:
                chains = compute_conserved_exon_chains(multi_exon_lists)
                for chain_idx, chain in enumerate(chains, start=1):
                    template = extract_sequence(fasta_path, chain.exons)
                    template_order_exons = _exons_in_template_order(chain.exons)
                    junctions = _compute_junction_positions(template_order_exons)
                    chain_id = f"{target_gene}_chain_{chain_idx}"
                    result.append(
                        ConservedExonChain(
                            id=chain_id,
                            exons=chain.exons,
                            template=template,
                            junction_positions_1based=junctions,
                            required_junction_positions_1based=junctions,
                        )
                    )
            except ValueError:
                # No conserved adjacencies — promote each transcript to its
                # own chain instead.
                for tl in multi_exon_lists:
                    template = extract_sequence(fasta_path, tl)
                    template_order_exons = _exons_in_template_order(tl)
                    junctions = _compute_junction_positions(template_order_exons)
                    tid = _tid_for_exons(tl)
                    result.append(
                        ConservedExonChain(
                            id=tid,
                            exons=tl,
                            template=template,
                            junction_positions_1based=junctions,
                            required_junction_positions_1based=junctions,
                        )
                    )

        # Add single-exon transcripts as junctionless chains
        for tl in single_exon_lists:
            template = extract_sequence(fasta_path, tl)
            tid = _tid_for_exons(tl)
            result.append(
                ConservedExonChain(
                    id=tid,
                    exons=tl,
                    template=template,
                    junction_positions_1based=[],
                    required_junction_positions_1based=[],
                )
            )

        return result


# ── helpers ──────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def _reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


_GTF_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def _parse_gtf_attributes(attr: str) -> dict[str, str]:
    """Parse the GTF attributes column into a key-value dict."""
    result: dict[str, str] = {}
    for match in _GTF_ATTR_RE.finditer(attr):
        key = match.group(1)
        value = match.group(2)
        result[key] = value
    return result
