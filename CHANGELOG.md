# Changelog

## v0.1.0 (2025-06-15)

Initial release of Brimer-PLAST, a local, open alternative to NCBI Primer-BLAST
that designs qRT-PCR primers spanning exon-exon junctions using primer3-py and
tnBLAST.

**Key features:**
- Primer design in two automatic modes (junction-spanning and intron-spanning)
- Specificity filtering via tnBLAST (genomic and transcriptome)
- Multi-target support (`--target-gene` and `--target-transcript`)
- PDF report with per-chain genome views
- Electron desktop GUI (beta)
- Cross-platform CI builds (macOS, Windows, Linux)