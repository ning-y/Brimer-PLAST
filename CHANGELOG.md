# Changelog

## v0.1.3 (2026-08-06)

Bug fix — primer pair names now report true amplicon coordinates.

**Changes:**
- Fix primer pair name so the end coordinate no longer counts the reverse
  primer twice (primer3's `PRIMER_RIGHT` is the 5'/outermost base of the
  reverse primer). Pair-name spans now match the reported product size.
- Extract a tested `make_pair_name()` helper; add `TestMakePairName`
  regression tests.
- Add bug ticket `docs/tickets/0001` and an archival correction note at the
  top of ADR 0002 (its recorded decision is left untouched).

## v0.1.2 (2026-07-05)

Electron production hardening — guard Python subprocess behind `app.isPackaged`
and fall back to `app.getVersion()` for the version string. CI and packaging
fixes for reliable tag-based version derivation.

**Changes:**
- Guard Python subprocess calls (`get-app-title`, `get-version`) behind
  `app.isPackaged` in production (Electron AppImage/package)
- Fall back to `app.getVersion()` when Python import fails in production
- Pass full version string through to electron-builder (PEP 440 → semver
  conversion with `.devN` → `-dev.N` fix)
- Use `fetch-depth: 0` instead of `fetch-tags` in CI checkout for correct
  setuptools-scm version derivation on tag pushes
- Add CI permission `contents: write` to release job for auto-release
- Update CONTEXT.md with versioning model, per-chain genome view, and
  contributing/non-contributing transcript concepts

## v0.1.1 (2026-06-15)

Versioning restructure — git tags are now the sole source of truth.

**Changes:**
- Reduce default `--num-return` from 50 to 10 for more focused results
- Auto-create GitHub Release with changelog on tag push
- Git tag becomes single source of truth for version; all hardcoded fallbacks
  replaced with sentinel `0.0.0`
- Electron installer version injected at build time via `build.mjs`,
  derived from Python `__version__`
- Update CONTEXT.md glossary to reflect versioning model

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