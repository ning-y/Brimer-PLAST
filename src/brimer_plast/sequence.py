"""Sequence extraction and coordinate mapping for Brimer-PLAST.

Extracts exon sequences from genome FASTA files, splices them into
transcript templates, and maps coordinates between template space and
genomic space.  Also provides the reverse_complement utility and
exon-ordering helper used throughout the codebase.
"""

from __future__ import annotations

from pathlib import Path

import pyfaidx

from brimer_plast.gtf import parse_gtf_all_transcripts
from brimer_plast.models import ExonInfo, GenomicFragment


# ── helper ───────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def exons_in_template_order(exons: list[ExonInfo]) -> list[ExonInfo]:
    """Return exons sorted in the order they appear in the spliced template."""
    if not exons:
        return []
    strand = exons[0].strand
    if strand == "-":
        return sorted(exons, key=lambda e: e.start, reverse=True)
    else:
        return sorted(exons, key=lambda e: e.start)


# ── sequence extraction ──────────────────────────────────────────────────────


def _extract_sequence_from_genome(
    genome: pyfaidx.Fasta,
    exons: list[ExonInfo],
) -> str:
    """Extract and concatenate exon sequences from an already-open FASTA handle.

    Like :func:`extract_sequence` but reuses an existing ``pyfaidx.Fasta``
    object instead of opening the file.  Designed for bulk transcriptome
    building where opening/closing per transcript is prohibitively slow.

    Raises:
        ValueError: If exons span different chromosomes, or seqid not found.
    """
    if not exons:
        return ""

    seqid = exons[0].seqid
    for exon in exons:
        if exon.seqid != seqid:
            raise ValueError(
                f"Exons span different chromosomes: {exon.seqid!r} != {seqid!r}. "
                "All exons for a target must share the same seqid."
            )

    if seqid not in genome:
        raise ValueError(
            f"Sequence {seqid!r} not found in FASTA genome. "
            f"Available sequences: {list(genome.keys())[:10]}"
        )

    ordered = exons_in_template_order(exons)

    parts: list[str] = []
    for exon in ordered:
        seq_record = genome[seqid]
        if seq_record is None:
            raise ValueError(f"Sequence {seqid!r} not found in FASTA genome.")
        sub = seq_record[exon.start - 1 : exon.end]
        if sub is None:
            raise ValueError(f"Subsequence {seqid}:{exon.start}-{exon.end} not found.")
        fragment = sub.seq
        if exon.strand == "-":
            fragment = reverse_complement(fragment)
        parts.append(fragment)

    return "".join(parts)


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
    try:
        return _extract_sequence_from_genome(genome, exons)
    finally:
        genome.close()


def build_transcriptome_fasta(
    fasta_path: str | Path,
    gtf_path: str | Path,
    output_path: str | Path,
) -> dict[str, list[ExonInfo]]:
    """Build a transcriptome FASTA by splicing all transcripts from GTF.

    Extracts exon sequences from the genome FASTA and splices them per
    transcript (reverse-complementing exons on the negative strand).
    Writes a FASTA file with one entry per transcript_id.

    Returns:
        The parsed transcript-to-exon map ``{transcript_id: [ExonInfo, ...]}``
        for all transcripts in the GTF.  Callers can use this to translate
        tnBLAST transcriptome coordinates back to genomic positions using
        the correct per-transcript exon list (see ``template_to_genomic``).
    """
    transcripts = parse_gtf_all_transcripts(gtf_path)
    genome = pyfaidx.Fasta(str(fasta_path), read_ahead=10_000)
    try:
        with open(output_path, "w") as f:
            for tid, exons in transcripts.items():
                template = _extract_sequence_from_genome(genome, exons)
                f.write(f">{tid}\n")
                # Write 60-char lines for readability
                for i in range(0, len(template), 60):
                    f.write(template[i : i + 60] + "\n")
    finally:
        genome.close()
    return transcripts


# ── coordinate mapping ───────────────────────────────────────────────────────


def template_to_genomic(
    template_pos_0based: int,
    primer_length: int,
    exons: list[ExonInfo],
) -> list[GenomicFragment]:
    """Map a template-relative position back to genomic coordinates.

    A single template position may map to multiple genomic fragments
    if it spans an exon-exon junction.

    Returns:
        List of :class:`GenomicFragment`.
    """
    ordered = exons_in_template_order(exons)
    fragments: list[GenomicFragment] = []

    current_template_pos = 0
    remaining_primer_len = primer_length
    primer_start_found = False

    for exon in ordered:
        exon_len = exon.end - exon.start + 1

        if not primer_start_found:
            if current_template_pos + exon_len > template_pos_0based:
                # Primer starts in this exon
                primer_start_found = True
                offset_in_exon = template_pos_0based - current_template_pos

                # Length available in this first exon
                len_in_this_exon = min(remaining_primer_len, exon_len - offset_in_exon)

                if exon.strand == "+":
                    g_start = exon.start + offset_in_exon
                    g_end = g_start + len_in_this_exon - 1
                else:
                    g_end = exon.end - offset_in_exon
                    g_start = g_end - len_in_this_exon + 1

                fragments.append(GenomicFragment(seqid=exon.seqid, start=int(g_start), end=int(g_end), strand=exon.strand))
                remaining_primer_len -= len_in_this_exon

            current_template_pos += exon_len
        else:
            # We already found the start, are there more bits in subsequent exons?
            if remaining_primer_len <= 0:
                break

            len_in_this_exon = min(remaining_primer_len, exon_len)

            if exon.strand == "+":
                g_start = exon.start
                g_end = g_start + len_in_this_exon - 1
            else:
                g_end = exon.end
                g_start = g_end - len_in_this_exon + 1

            fragments.append(GenomicFragment(seqid=exon.seqid, start=int(g_start), end=int(g_end), strand=exon.strand))
            remaining_primer_len -= len_in_this_exon

    return fragments


def genomic_range_to_fragments(
    g_start: int,
    g_end: int,
    exons: list[ExonInfo],
) -> list[GenomicFragment]:
    """Split a genomic coordinate range into exon-by-exon fragments.

    Intersects the range ``[g_start, g_end]`` (1-based inclusive) with each
    exon in *exons*, returning one :class:`GenomicFragment` per overlapping
    exon.  Exons are sorted in genomic order before processing.
    """
    ordered = sorted(exons, key=lambda e: e.start)
    fragments: list[GenomicFragment] = []
    for exon in ordered:
        o_start = max(g_start, exon.start)
        o_end = min(g_end, exon.end)
        if o_start <= o_end:
            fragments.append(GenomicFragment(seqid=exon.seqid, start=o_start, end=o_end, strand=exon.strand))
    return fragments
