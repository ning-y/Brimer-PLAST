"""GTF parsing for Brimer-PLAST.

Extracts exon coordinates and transcript/gene relationships from GTF
annotation files.  All coordinates are 1-indexed (GTF convention).

This module has no dependencies on other Brimer-PLAST modules except
``models``.  It does NOT open FASTA files or perform coordinate math.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from brimer_plast.models import ExonInfo

_GTF_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def _parse_gtf_attributes(attr: str) -> dict[str, str]:
    """Parse the GTF attributes column into a key-value dict."""
    result: dict[str, str] = {}
    for match in _GTF_ATTR_RE.finditer(attr):
        key = match.group(1)
        value = match.group(2)
        result[key] = value
    return result


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
                return attrs.get("gene_name") or attrs.get("gene") or attrs.get("gene_id")
    return None


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


# ── transcript → gene lookup ─────────────────────────────────────────────────


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
