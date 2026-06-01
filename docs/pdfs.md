# PDF output: context for tool search

## How dependencies work in this repo

Dependencies are declared in **three places**, depending on the deployment mode:

### 1. `pyproject.toml` (pip install)

The canonical Python dependency manifest. Only Python packages go here:

- `primer3-py >=2.1.0` — C-extension wrapping primer3
- `typer >=0.12.0,<2.0.0` — CLI framework
- `pyfaidx >=0.8.0` — indexed FASTA access

Optional dev deps: `pytest`, `ruff`.

Install via `pip install .` or `pip install -e ".[dev]"`.

**Constraints for adding a new dependency:**
- Must be a **Pure Python** package or a well-known C-extension (like reportlab, which ships its own C extensions).
- Must be installable from PyPI with a simple `pip install`.
- Noapt-get, no system libraries not already present in the Docker image (`python:3.12-slim` has `zlib1g` and `libgomp1`; the Dockerfile does not install additional system packages beyond these). If a new system library is required, the Dockerfile must be updated accordingly.

### 2. `flake.nix` (nix develop)

Defines the Nix dev shell. Currently builds:
- `tntblast` (custom derivation, minimal build)
- A Python environment (`python312.withPackages`) with `primer3`, `typer`, `pytest`, `pyfaidx`
- Plus `ruff`, `pip`

If a new Python dependency is added to `pyproject.toml`, the Nix shell must also be updated so `nix develop` provides it — either by adding it to the `withPackages` list, or by ensuring `pip install -e ".[dev]"` inside the venv (which the shellHook already does) can pull it from PyPI.

### 3. `Dockerfile` (docker build)

Two-stage OCI image:
- **Stage 1 (tntblast-builder):** `ubuntu:24.04` + `build-essential`, `g++`, `make`, `zlib1g-dev`, `wget` — compiles tnBLAST.
- **Stage 2 (final):** `python:3.12-slim` with `zlib1g` and `libgomp1` system packages. `pip install` runs on the project, and `tntblast` is copied from stage 1.

If a new system-level library is needed (e.g. libffi for cffi-based packages, libpango for WeasyPrint), it must be `apt-get install`ed in the final stage of the Dockerfile.

## Current dependencies (summary)

| Dependency | Version | Purpose | Type |
|---|---|---|---|
| python | >=3.12 | Runtime | Interpreter |
| primer3-py | >=2.1.0 | Primer design | PyPI |
| typer | >=0.12.0 | CLI | PyPI |
| pyfaidx | >=0.8.0 | FASTA indexing | PyPI |
| tntblast | 2.77 | Specificity filtering | Compiled C++ |
| pytest | >=8.0 | Testing | dev (PyPI) |
| ruff | >=0.5.0 | Linting | dev (PyPI/Nix) |

## Goal PDF description

Each invocation of Brimer-PLAST should produce a **PDF report** alongside (or instead of) the current terminal/TSV output. The report is archival — it goes in a lab notebook or supplementary material.

### Section 1 — Gene-view diagram

A visual representation of the target gene's structure with each filtered primer pair overlaid.

For each primer pair that passed specificity filtering, the diagram should show:

- Exons as **boxes** (proportional to length), connected by **lines** for introns
- The scale/ruler at the top showing genomic coordinates
- Forward primer → as an arrow above the exon/intron it binds to
- Reverse primer ← as an arrow (reverse-complement orientation)
- The amplicon as a bracket or shaded region connecting the two primers
- Label each pair (Pair 1, Pair 2, …) with its product size

If there are multiple conserved exon chains, each chain gets its own diagram.

### Section 2 — Primer pair table (same as current terminal output but as a PDF table)

| Pair | Forward (5→3) | Tm(°C) | %GC | Reverse (5→3) | Tm(°C) | %GC | Size (bp) |

### Section 3 — Experiment trace

A record of the invocation parameters that produced this report:
- Date and time
- Command-line arguments (genome, annotations, target, all primer design parameters)
- Timing/log lines from stderr (how many chains found, how many candidates designed, how many passed filtering)
- Version of Brimer-PLAST, primer3-py, tnBLAST
- Genome and annotation filenames (and checksums if feasible)

### Additional requirements

- **Pure Python library preferred**, or one that's easily available on PyPI with minimal system dependencies.
- Must render precise shapes (exon boxes, lines, arrows) with accurate proportional scaling.
- Must produce an archival-quality PDF (embedded fonts, vector graphics, no rasterisation of the diagram).
- Must work identically in: bare `pip install .`, `nix develop`, and the Docker image.

### What to search for

Candidate PDF libraries that meet the "pure Python" and "diagram drawing" needs. For each candidate, report:

1. **Library name and version**
2. **Install method** (pip package name, system deps if any)
3. **Capabilities relevant to the gene diagram** — can it draw rectangles, lines, arrows, text labels with exact positioning?
4. **Capabilities relevant to the table** — does it have a table/flowable framework for multi-page tables?
5. **Capabilities relevant to the experiment trace** — does it support multi-page text sections with headers?
6. **System library dependencies** (libffi, fontconfig, Pango, etc.) — what's needed beyond what's already in `python:3.12-slim`?
7. **Python version compatibility** (3.12+?)
8. **Example or prior art** — has it been used in bioinformatics reports before?
