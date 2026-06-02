"""Typed data models for Brimer-PLAST.

Provides dataclasses for the core data shapes used throughout the pipeline,
replacing generic ``dict[str, Any]`` with typed structures.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class GenomicFragment:
    """A contiguous region of the genome.

    Used to represent where a primer binds on the genome.  A primer that spans
    an exon-exon junction will have multiple GenomicFragments (one per exon).
    Coordinates are 1-based inclusive.
    """

    seqid: str
    start: int
    end: int
    strand: str  # "+" or "-"


@dataclass(frozen=True)
class ExonInfo:
    """Coordinates and strand of a single exon from a GTF annotation.

    Coordinates are 1-indexed (GTF convention).
    """

    seqid: str
    start: int
    end: int
    strand: str  # "+" or "-"


@dataclass
class PrimerPair:
    """A candidate primer pair designed by primer3-py.

    All numeric fields may be ``None`` if primer3 failed to compute them.

    *Primer3 template-relative positions*
        *forward_start*, *forward_len*, *reverse_start*, *reverse_len* are
        0-based start positions and lengths within the spliced template
        sequence, taken directly from primer3-py's output.

    *Genomic fragment lists*
        Both primer3 and tnBLAST-derived fragments are pre-computed in
        the CLI as ``list[GenomicFragment]``.  The renderer draws these
        directly without any coordinate math.
    """

    forward_seq: str = ""
    reverse_seq: str = ""
    forward_tm: float | None = None
    reverse_tm: float | None = None
    forward_gc: float | None = None
    reverse_gc: float | None = None
    product_size: int | None = None
    pair_penalty: float | None = None

    # primer3 template-relative positions (Panel A)
    forward_start: int | None = None  # 0-based in template
    forward_len: int | None = None
    reverse_start: int | None = None  # 0-based in template
    reverse_len: int | None = None

    # Pre-computed genomic fragments (primer3-derived, blue Panel A)
    primer3_forward_fragments: list[GenomicFragment] = field(default_factory=list)
    primer3_reverse_fragments: list[GenomicFragment] = field(default_factory=list)

    # Pre-computed genomic fragments (tnBLAST-derived, red Panel B)
    tnblast_forward_fragments: list[GenomicFragment] = field(default_factory=list)
    tnblast_reverse_fragments: list[GenomicFragment] = field(default_factory=list)

    # Which chain this pair belongs to
    chain_id: str = ""

    # Sequential number within the final filtered output (set by CLI)
    pair_number: int | None = None


@dataclass
class TargetInfo:
    """Deprecated. Replaced by :class:`ConservedExonChain`.

    Identification and spliced template sequence for a target gene/transcript.
    Kept for backward compatibility; the pipeline now returns
    ``list[ConservedExonChain]`` from :func:`get_target_information`.
    """

    id: str
    template: str
    exons: list[ExonInfo] = field(default_factory=list)


@dataclass
class ConservedExonChain:
    """A contiguous run of exons conserved across all transcripts of a gene.

    The *template* is the spliced mRNA (exons concatenated in genomic order).
    *junction_positions_1based* lists the 1-indexed positions of exon-exon
    boundaries within the template -- these are passed to Primer3 as
    ``SEQUENCE_OVERLAP_JUNCTION_LIST`` (soft penalty, not hard constraint).

    *required_junction_positions_1based* is the subset that at least one primer
    must actually overlap for the pair to survive post-filtering.  When empty,
    defaults to *junction_positions_1based* (all junctions are required).
    """

    id: str
    exons: list[ExonInfo]
    template: str
    junction_positions_1based: list[int] = field(default_factory=list)
    required_junction_positions_1based: list[int] = field(default_factory=list)


@dataclass
class GeneLocus:
    """The full genomic locus for a gene.

    Used for the Genome View in reports.
    """

    gene_name: str
    seqid: str
    strand: str
    transcripts: dict[str, list[ExonInfo]]  # tid -> exon list
    min_start: int
    max_end: int
