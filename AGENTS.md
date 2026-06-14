# Instructions for AI agents

This project is developed from scratch.

## Development principles

- **Test-driven development.** Write the test first, then implement. Every `src/` module should have a corresponding test file. Run `pytest` frequently.
- **Keep `CONTEXT.md` current.** When a new domain term enters the design, add it to the glossary.
- **Keep README.md current.** When CLI flags, output format, or behavior changes, update the user-facing documentation.

## Dependencies

Core project dependencies are provided by `flake.nix` at the project root. Use `nix develop --command <cmd>` to run commands inside the dev shell.

**Always run `pytest` and `python` through `nix develop --command`.** The `primer3-py` C-extension and `tntblast` binary are only available inside the nix shell. Do not invoke `.venv/bin/python` or `.venv/bin/pytest` directly — the `.venv/` directory is created by the nix `shellHook` and has `--system-site-packages` enabled; it is broken outside `nix develop`.

### Never skip integration tests

**Do not use `-k "not integration"` or any other filter that excludes integration tests.** The data is already present on disk (`tests/fixtures/ce11/`), so there is no reason to skip them. Doing so masks real failures that only manifest with real genome data.

Common commands:

```bash
# Run ALL tests (includes C. elegans integration, ~2 minutes)
nix develop --command timeout 240 pytest tests/

# Run a single test file (still runs if it happens to be an integration test)
nix develop --command pytest tests/test_pipeline.py -v

# Run a single test class or method
nix develop --command pytest tests/test_pipeline.py::TestPipelineResult -v

# Type-check with pyright (configured in pyproject.toml)
nix develop --command pyright

# Lint / format with ruff
nix develop --command ruff check src/
nix develop --command ruff format src/ --check
```

### Electron JS/HTML checks

**Before committing any change to `electron/`, run the JS syntax check:**

```bash
node --check electron/main.js
node --check electron/preload.js
node -e "
  const fs = require('fs');
  const html = fs.readFileSync('electron/renderer/index.html', 'utf8');
  const m = html.match(/<script>([\s\S]*?)<\/script>/);
  if (!m) { console.error('No <script> tag found in index.html'); process.exit(1); }
  new Function(m[1]);
  console.log('Render JS syntax OK');
"
```

This catches corrupted or invalid JavaScript before a bad file gets committed (the same type of corruption that `sed` can accidentally introduce).

Conda/mamba and its related tools are never allowed.

## Production packaging

The production distribution is an **OCI image** (Dockerfile).

- Build with: `docker build -t brimer-plast .`
- The image bundles Python, primer3-py, and a compiled tnBLAST.
- Users without Docker can install via `pip install .` if they have tnBLAST on PATH.

## Module organization

After a recent refactor, the monolithic `genome.py` and `cli.py` were split into focused modules:

```
src/brimer_plast/
├── __init__.py        # package exports (ConservedExonChain, ExonInfo, PrimerPair)
├── cli.py             # CLI entry point (typer), multi-target dispatch, --output-pdf
├── diagram.py         # PDF genome-view diagram drawing
├── filter.py          # tnBLAST result parsing and specificity filtering
├── genome.py          # re-exports from gtf.py + sequence.py; conserved chain detection
├── gtf.py             # GTF parsing (parse_gtf, build_transcript_to_gene_map)
├── log_config.py      # logging configuration
├── models.py          # data classes (ConservedExonChain, PrimerPair, ExonInfo, GeneLocus)
├── pdf_report.py      # PDF report generation (reportlab)
├── pipeline.py        # run_pipeline() — reusable core pipeline logic
├── primer.py          # primer3-py wrapper, default constants
├── sequence.py        # sequence extraction, coordinate conversion
└── tnblast.py         # tnBLAST subprocess wrapper, assay file writer
```

### CLI features added since initial design

- **Multi-target**: `--target-gene` and `--target-transcript` can be repeated:
  `--target-gene GAPDH --target-gene ACTB`. Each target runs independently.
- **PDF report**: `--output-pdf <path>` generates a PDF with genome views,
  per-chain diagrams, and filtered pair tables.
- **Pair naming**: Primer pairs are named `{short_tid}:{amplicon_start}-{amplicon_end}`
  (e.g. `9746.1:45-199`) instead of `pair_1`.
- **Fallback chains**: When no conserved exon-exon junctions exist across
  transcripts, the pipeline creates per-transcript chains (flagged with
  a user-visible warning).

## Integration tests

Integration tests use real C. elegans (WBcel235) genome data.

- Run `bash tests/fixtures/download-ce11.sh` to download (~42 MB).
- Tests are automatically skipped if the data isn't present.
- The download script is idempotent — running it multiple times just overwrites the files.

## Tool usage lessons

### When using the `edit` tool, match text exactly — do not extrapolate

The `edit` tool requires `oldText` to match the file byte-for-byte. Tabs and spaces are distinct. **Always read the exact target lines first** using the `read` tool with the correct `offset` and `limit`, then copy-paste that exact text into `oldText`. Do not construct `oldText` from memory, from `sed` output, or from `rg --context` — always use the `read` tool output.

Every edit failure in a long-running session came from guessing the indentation or reconstructing text rather than reading it verbatim first.

### After a series of edits, verify the file's state

Failed `edit` calls can still partially corrupt the file (e.g. appending garbage after a function's closing `}`). After a series of edits, run `git diff` or `read` the affected regions to confirm no corruption remains.

### Prefer `read` over `rg --context` for edit targets

`rg --context` shows matching lines but may display whitespace differently from the actual file. Always use `read` (which preserves the literal bytes) to get the text you paste into `edit`.

## Scopes

- **ci** — CI workflow definitions, build pipelines, artifact management
- **debug-archive** — Debug ZIP error-reporting system for Electron app
- **electron** — Electron main-process runtime, sidecar process lifecycle, and IPC message handling
- **electron-packaging** — Electron app bundling, electron-builder config, sidecar resource discovery
- **electron-ui** — Electron renderer UI, branding, header, theme, version display
