# Brimer-PLAST

Design qRT-PCR primers that span exon-exon junctions using primer3 and tnBLAST.

Named as a local, open alternative to [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/).

By default, every primer pair is constrained so that **at least one primer overlaps
an exon-exon junction** — this prevents amplification from genomic DNA contamination
in qRT-PCR experiments.  Use `--disable-junction-overlap` for genomic PCR.

## Quick start

```bash
# Docker
docker build -t brimer-plast .
docker run --rm -v /path/to/genomes:/data brimer-plast \
  --genome /data/genome.fna \
  --annotations /data/annotations.gtf \
  --target-gene GAPDH

# Native (requires tnBLAST on PATH and primer3-py installed)
pip install .
brimer-plast --genome genome.fna --annotations annotations.gtf --target-gene GAPDH

# Nix dev shell (includes tnBLAST + Python)
nix develop
brimer-plast --genome genome.fna --annotations annotations.gtf --target-gene homt-1
```

## Full usage

```
brimer-plast --genome <FASTA> --annotations <GTF> <--target-gene | --target-transcript> [options]
```

Required:
- `--genome` / `-g` — Genome FASTA file
- `--annotations` / `-a` — Gene annotation GTF file
- `--target-gene` or `--target-transcript` — Which gene/transcript to design primers for

Junction-spanning control:
- `--disable-junction-overlap` — Allow primers that do not span an exon-exon junction
  (use for genomic PCR rather than qRT-PCR)

Primer design options:
- `--num-return` / `-n` — Number of candidate pairs per conserved exon chain (default: 10)
- `--min-tm`, `--max-tm`, `--opt-tm` — Melting temperature range (default: 57/63/60 °C)
- `--min-size`, `--max-size`, `--opt-size` — Primer length (default: 18/25/20 bp)
- `--min-gc`, `--max-gc` — GC content percent (default: 40/60%)
- `--product-min`, `--product-max` — Amplicon size range (default: 100-400 bp)

Output options:
- `--tsv` — Tab-separated machine-readable output

## How it works

1. **Parse annotations** — Reads the GTF file, groups exons by `transcript_id`.
   For `--target-gene`, computes *conserved exon chains* (maximal contiguous runs
   of exons where every adjacent pair is present in all transcripts).
   For `--target-transcript`, uses that transcript's exons directly. It then
   determines which junctions are *unique* to this transcript (not shared
   with any sibling transcript of the same gene) for post-filtering.

2. **Extract template** — For each conserved exon chain, extracts and splices
   the exon sequences from the genome FASTA, producing a mature-mRNA template.

3. **Design primers** — Runs primer3 (via primer3-py C-extension, ~1000× faster
   than subprocess wrappers) on each template.  By default, passes the
   `SEQUENCE_OVERLAP_JUNCTION_LIST` parameter as a **soft penalty** to guide
   Primer3 toward junction-spanning pairs.  A **post-filter** then enforces the
   hard requirement: at least one primer (forward or reverse) must actually
   overlap a required junction.  For `--target-gene`, the required junctions
   are all conserved junctions.  For `--target-transcript`, the required
   junctions are only those unique to that transcript.  Candidates from all
   chains are pooled.

4. **Filter for specificity** — Runs tnBLAST to check each candidate pair for
   off-target amplification; only pairs with exactly one predicted amplicon
   (the on-target one) survive.

5. **Output** — Prints a ranked table of specificity-filtered primer pairs.

## When it errors

- **No conserved exon-exon junctions** — The target gene's transcripts share
  no common splice sites (or all transcripts have single exons).  Use
  `--disable-junction-overlap` to relax the constraint.
- **No unique junctions** (with `--target-transcript`) — The target transcript
  shares all its exon-exon junctions with other transcripts of the same gene,
  so no primer can be specific to this transcript alone.  Use
  `--disable-junction-overlap` to design non-junction-spanning primers (which
  may also amplify sibling transcripts).
- **No candidate primers** — Primer3 could not find qualifying pairs within
  the given temperature/size/product constraints.  Try widening the product
  size range (e.g. `--product-min 80 --product-max 300`).

## Limitations

- **Highly alternatively spliced genes** — When a gene has many transcripts
  that share no common exon-exon boundaries (e.g. CACNA1C with 34 isoforms),
  `--target-gene` reports no conserved exon chains.  Use
  `--target-transcript` with a specific transcript ID to target a single
  isoform.
- **Single-exon genes** — Genes with only one exon (including all
  mitochondrial genes) have no exon-exon junctions.  Use
  `--disable-junction-overlap` to design primers for these targets.

## Development

```bash
nix develop               # enter dev shell (Python 3.12 + compiled tnBLAST + ruff + pytest)
pytest tests/             # run all tests (93 tests: unit + integration)
bash tests/fixtures/download-ce11.sh  # download C. elegans for integration tests
```

## Output format

Human-readable (default):
```
Pair   Forward (5→3)                Tm(°C)   %GC   Reverse (5→3)                Tm(°C)   %GC   Size
----------------------------------------------------------------------------------------------------
1      TTCGTCGAAGGACTGCAGAC         60.0     55    TGCAGTGCTTTCGAGACCAT         60.0     50    281
```

TSV (`--tsv`):
```
pair	forward_seq	reverse_seq	forward_tm	reverse_tm	forward_gc	reverse_gc	product_size
1	TTCGTCGAAGGACTGCAGAC	TGCAGTGCTTTCGAGACCAT	60.0	60.0	55	50	281
```

## License

GPLv2 — same as primer3 and primer3-py.