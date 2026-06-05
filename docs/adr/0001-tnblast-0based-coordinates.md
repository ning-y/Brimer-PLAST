# tnBLAST amplicon ranges are 0-based inclusive

tnBLAST outputs `amplicon range = <start> .. <end>` where both start and end are **0-based inclusive** positions in the target database (e.g., genomic FASTA coordinates). Brimer-PLAST converts these to 1-based coordinates before storing them as ``GenomicFragment.start`` / ``GenomicFragment.end``.

## Context

`_compute_tnblast_fragments` in ``pipeline.py`` reads tnBLAST output coordinates and computes genomic fragments for the PDF diagram's Panel B (red). Getting this conversion wrong produces off-by-one primer placements that would look plausible but be biologically incorrect.

tnBLAST's output format is not formally documented, so we verified experimentally by constructing a known-sequence target, running tnBLAST with known primer positions, and reading the reported range.

## Experiment

FASTA database: a single 139-base sequence.
Assay: forward = bases 1-20 of the sequence, reverse = reverse complement of bases 120-139.

tnBLAST output::

    amplicon range = 0 .. 138
    amplicon length = 139

Forward primer is at 1-based positions 1-20 (0-based 0-19). The range `0 .. 138` covers the full sequence: 0-based inclusive 0 through 138 = 139 bases. A second test with the forward at 0-based position 50 gave `amplicon range = 50 .. 119`, confirming the convention.

## Conversion rule

| tnBLAST field | Meaning | Brimer-PLAST 1-based conversion |
|---|---|---|
| ``amplicon_start`` | 0-based start of forward primer binding site (positive strand) or 0-based start of reverse primer binding site (negative strand) | ``+ 1`` |
| ``amplicon_end`` | 0-based end of reverse primer binding site (positive strand) or 0-based end of forward primer binding site (negative strand) | ``+ 1`` |

The strand-reversal logic in `_compute_tnblast_fragments` (``pipeline.py`` lines 61-71) swaps forward/reverse fragment calculation when ``locus.strand == "-"`` because tnBLAST reports the amplicon in forward-strand genomic coordinates regardless of the gene's transcriptional strand.

## Status

Accepted.
