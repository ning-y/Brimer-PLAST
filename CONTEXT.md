# Brimer-PLAST

A command-line tool that designs qRT-PCR primer pairs for a target gene using
primer3, then filters out non-specific primers using tnBLAST. Named as a local,
open alternative to NCBI's Primer-BLAST.

## Language

**Brimer-PLAST**:
The tool itself. A wrapper that orchestrates primer3 (via primer3-py) and tnBLAST
to produce specificity-filtered primer pairs from a user-supplied genome, annotation,
and target gene.
_Avoid_: Primer-BLAST, the tool, the program

**qRT-PCR (quantitative reverse-transcription PCR)**:
The intended use case for Brimer-PLAST. Primers are designed to span exon-exon
junctions so that they amplify only spliced mRNA, not contaminating genomic DNA.
_Avoid_: RT-PCR, real-time PCR

**Exon-exon junction**:
The boundary between two adjacent exons in a spliced mRNA. Brimer-PLAST requires
at least one primer in each pair to overlap this boundary (default behavior), which
prevents amplification from genomic DNA.
_Avoid_: Splice junction, exon boundary

**Junction-spanning primer**:
A primer whose sequence straddles an exon-exon junction — the primer has one
portion in the upstream exon and the other in the downstream exon. This is the
default constraint for all primers Brimer-PLAST designs.
_Avoid_: Junction-overlapping primer, spanning primer

**--disable-junction-overlap**:
The CLI flag to turn off the junction-spanning requirement. Use for genomic PCR
where gDNA contamination is not a concern.
_Avoid_: --genomic-pcr, --no-junction

**Target gene**:
A gene name (e.g. GAPDH) specified with ``--target-gene``. Brimer-PLAST groups
its exons by transcript ID and identifies conserved exon-exon junctions present
in all transcripts of that gene.

Multiple targets can be specified by repeating the flag:
``--target-gene GAPDH --target-gene ACTB``. Each target is processed as an
independent invocation, producing its own primer set and report.
_Avoid_: Gene name, gene symbol

**Target transcript**:
A specific transcript ID (e.g. NM_001289746.1) specified with
``--target-transcript``. Brimer-PLAST uses all splice junctions of that one
transcript for the overlap requirement.

Multiple targets can be specified by repeating the flag:
``--target-transcript NM_001289746.1 --target-transcript NM_001234567.1``.
Each target is processed as an independent invocation.
_Avoid_: Transcript ID, refseq ID

**Conserved exon chain**:
A maximal contiguous run of exons where every adjacent pair's shared boundary
appears in all transcripts of the gene. Each conserved chain becomes its own
spliced template for primer design. The report includes a full sequence view
of each chain, organized by exon in 5' to 3' biological order.
_Avoid_: Exon block, conserved run

**Genomic coordinate notation**:
Brimer-PLAST uses 1-based coordinates. In report headers and technical views,
exons on the negative strand are notated as `high-low` (e.g., `chr1:2000-1000`)
to reflect the 5' to 3' orientation of the transcript, while positive strand 
exons use standard `low-high` (e.g., `chr1:1000-2000`).
_Avoid_: 0-based coordinates, BED format

**Candidate primer pair**:
A forward and reverse primer produced by primer3, prior to specificity filtering
by tnBLAST. These may amplify off-target regions of the genome.
_Avoid_: Raw primers, unfiltered primers

**Forward primer**:
The primer that is identical in sequence to a segment of the mRNA (sense strand). 
In the spliced template, it is always 5' of the reverse primer. Genomically, for 
positive-strand genes, it has smaller coordinates than the reverse primer; for 
negative-strand genes, it has larger coordinates.
_Avoid_: Left primer, sense primer

**Reverse primer**:
The primer that is the reverse-complement of a segment of the mRNA. It binds to 
the cDNA (antisense strand) during PCR. In the spliced template, it is always 
3' of the forward primer.
_Avoid_: Right primer, antisense primer

**Specificity-filtered primer pair**:
A candidate primer pair that survived tnBLAST's screen — i.e., no off-target
amplicons were predicted at the given annealing temperature.
_Avoid_: Good primers, clean primers

**Off-target amplification**:
An amplicon predicted by tnBLAST at a locus other than the intended target gene.
The presence of off-target amplification disqualifies a candidate primer pair.
_Avoid_: Non-specific binding, cross-reactivity

**Primer set**:
The complete collection of specificity-filtered primer pairs for a given target gene
that Brimer-PLAST returns to the user.
_Avoid_: Primer pool, primer collection

**primer3-py**:
A Python C-extension library that wraps primer3's C code directly. Brimer-PLAST
calls it via ``import primer3`` rather than shelling out to ``primer3_core``.
~1000× faster than subprocess wrappers.
_Avoid_: primer3_core, primer3 subprocess

**SEQUENCE_OVERLAP_JUNCTION_LIST**:
A Primer3 sequence tag that specifies 1-based positions where at least one primer
(forward or reverse) must overlap the junction by at least
``PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION`` nucleotides.
_Avoid_: SEQUENCE_TARGET, overlap positions, junction list

**tnBLAST**:
ThermonucleotideBLAST — a C++ tool that searches DNA databases with
physically relevant measures (free energy, melting temperature) rather than
BLAST's heuristic similarity scores. Brimer-PLAST builds it in minimal mode
(FASTA input, no MPI/NCBI toolkit).
_Avoid_: BLAST, NCBI BLAST

**Transcript-specificity**:
When ``--target-transcript`` is specified, Brimer-PLAST compares the target
transcript's exon-exon junctions against all other transcripts of the same
sibling gene. A junction is *unique* only if its adjacency pair does not
appear in any sibling transcript. Primer3 receives all junctions as a soft
penalty ("SEQUENCE_OVERLAP_JUNCTION_LIST"), then a post-filter drops pairs
where neither primer spans a unique junction.
_Avoid_: Isoform-specific, splice-variant-specific

**Unique junction**:
An exon-exon adjacency pair that appears in the target transcript but not in
any other transcript of the same gene. Used as the post-filter criterion for
``--target-transcript`` mode.
_Avoid_: Private junction, distinct junction

**Required junctions**:
The minimum set of 1-based template positions where at least one primer in
each pair must actually overlap the junction. For ``--target-gene`` this
equals *all conserved junctions*. For ``--target-transcript`` this equals
only the *unique junctions*. Stored as
``ConservedExonChain.required_junction_positions_1based``.
_Avoid_: Filtered junctions, overlap targets

**Post-filter (junction overlap)**:
After Primer3 returns candidate pairs, ``design_primers`` checks each pair's
forward and reverse primer coordinates against the required junction
positions. A primer at 0-based start *s* with length *L* spans a 1-based
junction *j* iff ``s < j - 1 < s + L``. Pairs where neither primer spans any
required junction are dropped. This is a hard enforcement on top of
Primer3's soft penalty.
_Avoid_: Hard constraint, secondary filter

**ConservedExonChain model**:
A typed data structure (``ConservedExonChain`` dataclass) that holds a chain ID,
the list of constituent exons, the spliced template string, 1-based junction
positions (Primer3), required junction positions (post-filtering), and naming
metadata (``transcript_offset``, ``representative_tid``, ``short_tid_length``).
A ``fallback`` flag marks chains created from a single transcript when no
conserved exon-exon junctions exist across all transcripts of the gene.
_Avoid_: TargetInfo, chain model

**Primer pair model**:
A typed data structure (``PrimerPair`` dataclass) that holds forward/reverse
sequences, Tm, GC%, product size, penalty, a descriptive ``pair_name``
(``{short_tid}:{amplicon_start}-{amplicon_end}``), and pre-computed genomic
fragment lists for both primer3 (Panel A) and tnBLAST (Panel B) sources.
_Avoid_: Raw dict, untyped pair

**Genomic fragment**:
A typed data structure (``GenomicFragment`` dataclass) describing a single
contiguous region of the genome with seqid, start, end, and strand. A
junction-spanning primer will have multiple fragments — one per exon.
Fragments are pre-computed during pipeline execution and stored on
``PrimerPair`` so that the renderer never performs coordinate math.
_Avoid_: Tuple, raw coordinate, (seqid, start, end, strand)

**Exon info model**:
A typed data structure (``ExonInfo`` dataclass) that holds seqid, start, end,
and strand for a single exon parsed from GTF annotations.
_Avoid_: Raw dict, untyped exon coordinate

**Product size range**:
The acceptable PCR amplicon length (bp) passed to primer3 as
``PRIMER_PRODUCT_SIZE_RANGE``. Controlled via ``--product-min`` and
``--product-max``. Defaults to 80–200 bp for qRT-PCR.
_Avoid_: Amplicon size range, fragment length

**Fallback chain**:
When no conserved exon-exon junctions exist across all transcripts of a gene,
Brimer-PLAST creates one ``ConservedExonChain`` per transcript as a fallback
(``ConservedExonChain.fallback == True``). Primers from fallback chains may
not amplify all transcripts of the gene. A warning is printed to stderr.
_Avoid_: Degenerate chain, single-transcript chain, per-transcript chain

**Gene Locus**:
The full genomic coordinates and collection of all annotated transcripts for a 
given gene. Used as the biological context for report visualizations.

**Genome View**:
A visualization in the PDF report showing the gene locus, coordinate axis, all 
stacked transcripts, and designed primer binding sites. The coordinate axis 
always increases from left to right (genomic order), meaning negative-strand 
genes will appear to be oriented "backwards" (5' on the right). Replaces the 
internal design-chain visualization.
_Avoid_: Chromosome flip, reversed axis

**Transcriptome-to-genome mapping**:
The process of translating coordinates from a spliced mRNA sequence back to 
genomic chromosome coordinates. Used to ensure compatibility between 
different data sources (Primer3, tnBLAST) in the Genome View.

**Primer pair name**:
A descriptive identifier for each primer pair in the format
``{short_tid}:{amplicon_start}-{amplicon_end}``, where *amplicon_start*
and *amplicon_end* are 1-based positions relative to the mature mRNA
sequence of the representative transcript for that chain.
_Avoid_: pair_1, pair_2, generic counter

**Short transcript ID (short_tid)**:
The shortest suffix of a transcript ID (plus 4 characters) that uniquely
identifies every transcript in the gene. Computed per-gene: find L such
that the last L characters of all transcript IDs are unique, then use
L+4 (the fifth-shortest suffix). Canonical example:
``NM_001289746.1`` with a sibling ``NM_001234567.1`` has L=3 (``6.1``
vs ``7.1``), so short_tid_length=7 → ``9746.1``.
_Avoid_: short_name, abbreviation, short form

**Representative transcript**:
For each conserved exon chain, the alphanumerically first transcript
ID (lexicographic order) among those whose exons contain the chain's
exons. Used to compute the chain's transcript offset and to derive
short_tid for the primer pair name.
_Avoid_: Canonical transcript, default transcript

**Transcript offset (chain_offset)**:
The cumulative length (in bases) of the representative transcript's
exons that appear 5' of the conserved chain's first exon, measured in
mRNA 5'→3' order (TSS = position 1). Added to template-relative
amplicon coordinates to produce transcript-relative coordinates.
Stored on ``ConservedExonChain.transcript_offset``.

**Contributing transcript**:
For a given conserved exon chain, a transcript whose exon set (by genomic
start/end coordinates) contains every exon in the chain. In that chain's
genome view, contributing transcripts are drawn with normal fill and solid
stroke.
_Avoid_: Active transcript, relevant transcript

**Non-contributing transcript**:
For a given conserved exon chain, a transcript whose exon set does NOT
contain every exon in the chain. In that chain's genome view,
non-contributing transcripts are drawn with no fill, dashed exon borders,
and dashed intron arrows to clearly distinguish them from the contributing
transcripts.
_Avoid_: Irrelevant transcript, inactive transcript

**Per-chain genome view**:
One genome view per conserved exon chain in the PDF report, rather than one
view for the whole target. Each per-chain view shows all transcripts of the
gene (colored by contribution status) and only the primer pairs designed
from that specific chain. Chains with zero filtered pairs still show an
explicit message. Per-chain views are concatenated in flat order (Option A).
_Avoid_: Single-view report, combined view

## Example dialogue

**Dev**: I ran Brimer-PLAST on the mouse genome targeting the GAPDH gene. It returned a primer set of 12 specificity-filtered primer pairs.

**Domain expert**: So these 12 pairs all passed tnBLAST's screen — no off-target amplification predicted at the default annealing temperature?

**Dev**: That's right. The candidate primer pairs that did show off-target hits were dropped. I'm left with the specificity-filtered ones.

**Domain expert**: Good — that's what we'd call working primers. You can order synthesis for the top 3.

**Dev**: I had to use --disable-junction-overlap because GAPDH has only one conserved exon chain with two exons, and primer3 couldn't find any junction-spanning candidates within the product size range.

**Domain expert**: That happens when exons are too short. Try widening the product size range or reduce the number of junctions required.