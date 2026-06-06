# Primer pair names use fifth-shortest disambiguating transcript suffix

Primer pairs are named ``{short_tid}:{amplicon_start}-{amplicon_end}``, where
*short_tid* is the fifth-shortest suffix that uniquely identifies every
transcript in the gene, and coordinates are 1-based relative to the mature
mRNA of the representative transcript for that conserved exon chain.

## Context

Brimer-PLAST originally named primer pairs ``pair_1``, ``pair_2``, etc. — a
flat counter with no biological meaning.  Users requested names that allow
them to identify which transcript a pair was designed on and where the
amplicon lies, so they can cross-reference with lab notebooks, PCR plates,
and external databases without opening the PDF report each time.

Several naming schemes were considered before settling on the current one.

## Considered alternatives

### A. Last-4-digits of accession number

Extract the last four digits of the numeric portion of known-accession
transcripts (RefSeq, Ensembl), appending a version suffix if present.
E.g., ``NM_001289746.1`` → ``9746.1``, ``NM_001289746`` → ``9746``.

*Rejected because*: it ties the naming convention to a specific set of
annotation sources (RefSeq, Ensembl).  Users running the tool on custom
genomes (FlyBase, WormBase, SGD, TAIR) would get opaque fallback
behaviour.  The third-shortest-suffix approach works for any set of
transcript IDs with no source-specific patterns.

### B. Full transcript ID

Use the complete transcript name (e.g., ``NM_001289746.1:45-199``).

*Rejected because*: transcript IDs are long (20+ characters), making names
unwieldy in lab notebooks, spreadsheets, and column headers.  A suffix
is sufficient to disambiguate all transcripts in a single gene (which
typically has 2–20 isoforms), and adding two extra characters accounts
for the user's mental model of "just the distinguishing bit."

### C. Shortest possible suffix

Use the smallest L such that last-L characters of all transcript IDs
in the gene are unique.  E.g., if L=3 already disambiguates, use that.

*Rejected because*: a 3-character suffix (``6.1``, ``7.2``) feels cryptic
and error-prone for manual transcription.  The fifth-shortest (L+4) gives
enough context to be recognisable without being verbose.

## Decision

For each target gene, compute the shortest suffix length L such that the
last L characters of every transcript ID are unique within that gene.
Then use ``short_tid_length = L + 4``.  This is "fifth-shortest" because
L is the first shortest, L+1 is the second, L+2 is the third, L+3 is the
fourth, and L+4 is the fifth.

For each conserved exon chain, pick the alphanumerically first transcript
(lexicographic order) among those whose exon keys appear in the chain.
This transcript supplies the short_tid and the exon-to-TSS offset for
coordinate translation.

The pair name is then:
``{short_tid}:{chain_offset + forward_start + 1}-{chain_offset + reverse_start + reverse_len}``

All coordinates are 1-based relative to the mature mRNA TSS (position 1),
regardless of the gene's chromosomal strand.

## Consequences

Positive:
- Names are deterministic, stable across runs, and bi-directionally
  traceable back to the annotation.
- No dependency on accession-pattern matching — works with any genome
  annotation.
- The L+4 heuristic produces recognisable short IDs (typically 7–9
  characters) that fit in spreadsheets and column headers.

Negative:
- Adding four characters when L=1 already works feels gratuitous but is
  intentional (see "Considered alternatives C").
- If the annotation is later updated with new transcript IDs that happen
  to share a long suffix with existing ones, the short_tid_length may
  increase, changing all existing pair names.

## Status

Accepted.