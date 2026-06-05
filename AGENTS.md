# Instructions for AI agents

This project is developed from scratch.

## Development principles

- **Test-driven development.** Write the test first, then implement. Every `src/` module should have a corresponding test file. Run `pytest` frequently.
- **Phase order is strict.** Do not skip ahead. Each phase depends on the previous one. If a phase isn't checked off in TODO.md, don't start it.
- **Keep `CONTEXT.md` current.** When a new domain term enters the design, add it to the glossary.

## Dependencies

Core project dependencies are provided by `flake.nix` at the project root. Use `nix develop --command <cmd>` to run commands inside the dev shell.

**Always run `pytest` and `python` through `nix develop --command`.** The `primer3-py` C-extension and `tntblast` binary are only available inside the nix shell. Do not invoke `.venv/bin/python` or `.venv/bin/pytest` directly — the `.venv/` directory is created by the nix `shellHook` and has `--system-site-packages` enabled; it is broken outside `nix develop`.

Common commands:

```bash
# Run all tests (includes C. elegans integration, ~2 minutes)
nix develop --command timeout 240 pytest tests/

# Run a single test file
nix develop --command pytest tests/test_pipeline.py -v

# Run a single test class or method
nix develop --command pytest tests/test_pipeline.py::TestPipelineResult -v
```

Conda/mamba and its related tools are never allowed.

## Production packaging

The production distribution is an **OCI image** (Dockerfile).

- Build with: `docker build -t brimer-plast .`
- The image bundles Python, primer3-py, and a compiled tnBLAST.
- Users without Docker can install via `pip install .` if they have tnBLAST on PATH.

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
