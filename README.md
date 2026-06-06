# Brimer-PLAST

Design qRT-PCR primers that span exon-exon junctions using primer3 and tnBLAST.

Named as a local, open alternative to [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/).

Primers are designed in two automatic modes:
- **Junction mode** — at least one primer overlaps an exon-exon junction
- **Intron mode** — forward and reverse primers fall on different exons with
  a total intronic separation >1000 bp
Results are deduplicated by primer sequence, with junction pairs given priority.

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
  (repeatable: `--target-gene GAPDH --target-gene ACTB` for multiple targets)

Primer design options:
- `--num-return` / `-n` — Number of candidate pairs per conserved exon chain (default: 50)
- `--min-tm`, `--max-tm`, `--opt-tm` — Melting temperature range (default: 57/63/60 °C)
- `--min-size`, `--max-size`, `--opt-size` — Primer length (default: 18/25/20 bp)
- `--min-gc`, `--max-gc` — GC content percent (default: 40/60%)
- `--product-min`, `--product-max` — Amplicon size range (default: 80-200 bp)
- `--max-amplicon` — Maximum tnBLAST amplicon search length (default: 2000)

Output options:
- `--output-pdf <path>` — Write a PDF report with per-chain genome views.
  Repeat once per target for multi-target runs
  (e.g. `--output-pdf gapdh.pdf --output-pdf actb.pdf`).
  If omitted, a report is auto-generated; use `--no-pdf` to suppress.
- `--no-pdf` — Suppress PDF report generation entirely.
- `--verbose` / `-v` — Increase verbosity. `-v` for pipeline progress,
  `-vv` for per-pair fragment-list details and template coordinates.
- `--tsv` — Tab-separated machine-readable output

## How it works

1. **Parse annotations** — Reads the GTF file, groups exons by `transcript_id`.
   For `--target-gene`, computes *conserved exon chains* (maximal contiguous runs
   of exons where every adjacent pair is present in all transcripts).
   If no conserved junctions exist across the gene, creates one chain per
   transcript as a **fallback** (with a user-visible warning).
   For `--target-transcript`, uses that transcript's exons directly. It then
   determines which junctions are *unique* to this transcript (not shared
   with any sibling transcript of the same gene) for post-filtering.

2. **Extract template** — For each conserved exon chain, extracts and splices
   the exon sequences from the genome FASTA, producing a mature-mRNA template.

3. **Design primers** — Runs primer3 (via primer3-py C-extension, ~1000× faster
   than subprocess wrappers) on each template in two modes:

   - **Junction mode** — passes `SEQUENCE_OVERLAP_JUNCTION_LIST` as a soft
     penalty to guide Primer3 toward junction-spanning pairs, then hard
     post-filters: at least one primer must overlap a required junction.
   - **Intron mode** — no junction constraint; post-filters so forward and
     reverse primers are on different exons with total intronic separation
     >1000 bp.

   For `--target-gene`, the required junctions are all conserved junctions.
   For `--target-transcript`, the required junctions are only those unique
   to that transcript.  Each mode requests `num_return + 500` from Primer3,
   ensuring enough candidates for the downstream merge.  Results are
   deduplicated by primer sequence with junction pairs given priority, then
   truncated at the user's `--num-return` cap.  Each pair is assigned a
   descriptive name (`{short_tid}:{amplicon_start}-{amplicon_end}`, e.g.
   `9746.1:45-199`).

4. **Filter for specificity** — Runs tnBLAST to check each candidate pair for
   off-target amplification; only pairs with exactly one predicted amplicon
   (the on-target one) survive.

5. **Output** — Prints a ranked table of specificity-filtered primer pairs.
   Optionally writes a PDF report with per-chain genome views, transcript
   contribution styling, and filtered pair tables.

## When it errors

- **No candidate primers** — Neither junction nor intron mode could find
  qualifying pairs within the given constraints.  This can happen with
  single-exon targets (no junctions, no introns) or when constraints are
  too tight.  Try widening the product size range (e.g. `--product-min 80
  --product-max 300`).
- **No unique junctions** (with `--target-transcript`) — The target transcript
  shares all its exon-exon junctions with other transcripts of the same gene,
  so no junction-spanning primer can be specific to this transcript alone.
  The intron mode may still produce pairs, but they won't be isoform-specific.

## Limitations

- **Highly alternatively spliced genes** — When a gene has many transcripts
  that share no common exon-exon boundaries (e.g. CACNA1C with 34 isoforms),
  `--target-gene` reports no conserved exon chains.  Use
  `--target-transcript` with a specific transcript ID to target a single
  isoform.
- **Single-exon genes** — Genes with only one exon (including all
  mitochondrial genes) have no junctions and no introns.  Neither design
  mode can produce primers for these targets.

## Development

```bash
nix develop               # enter dev shell (Python 3.12 + compiled tnBLAST + ruff + pytest + pyright)
pytest tests/             # run all tests (190 tests: unit + integration)
pyright                   # type-check the codebase
ruff check src/           # lint
bash tests/fixtures/download-ce11.sh  # download C. elegans for integration tests
```

## Output format

Human-readable (default):
```
Pair Name            Forward (5'→3')             Tm(°C)    %GC   Reverse (5'→3')             Tm(°C)    %GC   Size
-------------------------------------------------------------------------------------------------------------------
9746.1:45-199        TTCGTCGAAGGACTGCAGAC         60.0      55    TGCAGTGCTTTCGAGACCAT         60.0      50    281
                                                         (TGCAGTGCTTTCGAGACCAT)
```

The first line shows the reverse primer as stored (same strand as forward).
The indented second line shows the reverse-complement form (actual PCR binding strand).

TSV (`--tsv`):
```
pair_name	forward_seq	reverse_seq	reverse_rc	forward_tm	reverse_tm	forward_gc	reverse_gc	product_size
9746.1:45-199	TTCGTCGAAGGACTGCAGAC	TGCAGTGCTTTCGAGACCAT	TGCAGTGCTTTCGAGACCAT	60.0	60.0	55	50	281
```

## License

GPLv2 — same as primer3 and primer3-py.