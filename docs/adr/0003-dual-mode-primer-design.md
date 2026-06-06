# 0003 — Dual-mode primer design (junction + intron) replaces —disable-junction-overlap

The `--disable-junction-overlap` flag was removed. Instead, the tool always
runs two Primer3 passes on the same spliced template: one requiring at least
one primer to span an exon-exon junction (mode A), and one requiring the
forward and reverse primers to fall on different exons with a total intronic
separation >1000 bp (mode B). Results are deduplicated by primer sequence,
prioritising mode A over mode B, and truncated at the user's `--num-return` cap.

## Considered options

- **Single flag (`--disable-junction-overlap`)** — forced the user to choose
  between strategies up front; missed cases where one strategy failed but the
  other would work.
- **Genomic FASTA template for mode B** — would require building a second
  template with introns included; more complex coordinate mapping. Rejected
  in favour of using the same spliced template and post-filtering on
  genomic intron distances computed from the exon list.
- **Single-call, post-filter only** — run Primer3 once with no constraint,
  then split by which filter each pair passes. Risks Primer3 spending its
  budget on pairs that fail both filters.

## Consequences

- Single-exon targets (snoRNAs, miRNAs, single-exon protein-coding genes)
  produce zero pairs from both modes and get a warning + empty output.
- The `--disable-junction-overlap` flag no longer exists — users who relied
  on it for single-exon work must accept the empty result.
- `PRIMER_NUM_RETURN` passed to Primer3 is now `num_return + 500` per mode
  to ensure the dedup step has enough candidates to fill the final cap.