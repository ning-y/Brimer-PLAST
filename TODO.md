# TODO: Electron desktop app for Brimer-PLAST

Turn the command-line tool into a double-clickable desktop app (`.dmg` for Mac,
`.exe` for Windows) that biologists can run offline with no internet, no server,
and no authentication.

**Approach**: monorepo. The Electron wrapper lives alongside the Python library.

---

## Final repo layout

```
brimer-plast/
├── electron/                        # NEW — Electron app
│   ├── main.js                      #   Main process: spawns sidecar, IPC bridge
│   ├── preload.js                   #   Safe context bridge for renderer
│   ├── renderer/
│   │   ├── index.html               #   GUI: full HTML form + results table
│   │   └── style.css                #   (optional — inline is fine too)
│   ├── package.json                 #   npm: electron + electron-builder
│   └── electron-builder.yml         #   Build config for .dmg / .exe
├── sidecar.py                       # NEW — stdin/stdout JSON-RPC wrapper
├── pyproject.toml                   # EXISTING (add pyinstaller to dev deps)
├── .github/workflows/
│   ├── ci.yml                       # EXISTING — pytest, pyright, ruff
│   └── desktop.yml                  # NEW — build tnBLAST + PyInstaller + Electron
├── src/brimer_plast/                # EXISTING — Python library
├── tests/                           # EXISTING
├── Dockerfile                       # EXISTING
├── flake.nix                        # EXISTING
├── README.md                        # EXISTING
├── CONTEXT.md                       # EXISTING
├── AGENTS.md                        # EXISTING
└── TODO.md                          # THIS FILE
```

---

## Phase 0 — Cross-compile dependencies

Goal: produce a `tntblast` binary for each target platform. No submodule, no
fork — CI downloads the upstream source, compiles it, and discards it.

### 0.1 — tnBLAST for macOS (Intel x86_64 + ARM64)

Both architectures use the same commands — just different GitHub runners
(`macos-13` for Intel, `macos-14` for ARM).

```bash
# On a macOS runner:
curl -sL -o tnblast.tar.gz \
  https://github.com/jgans/thermonucleotideBLAST/archive/refs/tags/v2.77.tar.gz
tar xzf tnblast.tar.gz --strip-components=1 -C tnblast-src
cd tnblast-src

# Minimal build: no MPI, no NCBI toolkit, no OpenMP
make \
  CC=clang++ \
  FLAGS="-O3 -Wall -std=c++14" \
  LIBS="-lm -lz" \
  -j$(sysctl -n hw.logicalcpu)
```

- `CC=clang++` — macOS runners ship Apple Clang, not GCC. Clang accepts the
  same flags for this code.
- No `-fopenmp` — Apple Clang doesn't ship OpenMP without libomp. Single-
  threaded tnBLAST is fine for our workload.
- `options.cpp` uses `getopt_long()` — macOS has it in `<getopt.h>` (BSD
  systems). ✅
- All MPI code is guarded by `#ifdef USE_MPI` — removed by omission. ✅
- zlib is available on macOS by default. ✅

**Unknown**: None. This should compile cleanly. If it doesn't, the only likely
cause is a missing `-lz` (zlib not in the default linker path).

**Verification**: `file tntblast` should show the correct architecture
(x86_64 or arm64). `./tntblast --help` should print usage.

### 0.2 — tnBLAST for Windows (x86_64)

Windows runners don't have GCC. We install MSYS2 + MinGW-w64 via a GitHub
Action, then build inside it.

```bash
# On a windows-latest runner, after msys2/setup-msys2@v2:
pacman -S --noconfirm mingw-w64-x86_64-gcc mingw-w64-x86_64-zlib
cd tnblast-src
mingw32-make \
  CC=g++ \
  FLAGS="-O3 -Wall -std=c++14" \
  LIBS="-lm -lz"
```

- `mingw32-make` is MinGW's `make`.
- MinGW provides `getopt.h` and zlib. ✅
- Output is `tntblast.exe`.

**Unknown**: Whether MinGW's `getopt.h` fully implements `getopt_long()`.
tnBLAST uses it with the standard `option` struct and `long_index`. MinGW's
implementation should handle this, but it's the single point where Windows
could differ.

**Fallback if MinGW getopt is incomplete**: Bundle a standalone `getopt.h`
implementation (there are public-domain ones in ~50 lines). Add it to the
CI build script as a header override.

### 0.3 — primer3-py from source on macOS ARM

Prebuilt wheels on PyPI cover:
- ✅ Linux x86_64 (cp312)
- ✅ macOS Intel x86_64 (cp312)
- ✅ Windows amd64 (cp312)
- ⬜ macOS ARM64 (cp312) — only cp39 has an ARM64 wheel

On ARM runners (`macos-14`), `pip install primer3-py` will try the prebuilt
wheel, find no cp312 ARM match, and fall back to building from source. The
bundled `libprimer3` C code is standard C89; the Cython extensions compile
with Apple Clang. This **should** work, but it's the highest-risk item in
Phase 0.

```bash
pip install primer3-py  # builds from source on ARM
```
Then:
```bash
python -c "import primer3; primer3.design_primers(...)"  # smoke test
```

**Unknown**: Does `libprimer3/Makefile` (which uses `gcc`) work with Apple
Clang on ARM? On macOS, `gcc` is symlinked to `clang`. The Makefile uses
`-Wno-deprecated`, `-O2`, `-Wall`, `-D__USE_FIXED_PROTOTYPES__` — all
standard flags that Clang accepts. But the `Makefile.OSX` file exists for a
reason (it uses `libtool` instead of `ar`, sets SDKROOT, etc.). The main
Makefile may or may not produce correct static libraries on ARM macOS.

**If it fails**: The `setup.py` already handles Windows with special flags
(`TESTOPTS=--windows`). We may need to add an equivalent macOS-ARM path.

---

## Phase 1 — Python sidecar

The sidecar is a thin bridge. Electron spawns it as a child process and talks
to it over stdin/stdout JSON. No network, no temp files.

### 1.1 — Write `sidecar.py` (~100 lines)

Location: repo root (`/home/ning/repos/Brimer-PLAST/sidecar.py`)

Protocol — JSON lines, one request per line on stdin, one response per line
on stdout. Request/response pairs matched by `id` field:

**Request** (one JSON line from Electron → sidecar):
```json
{
  "id": 1,
  "command": "run_pipeline",
  "params": {
    "genome": "/Users/alice/genomes/hg38.fna",
    "annotations": "/Users/alice/genomes/hg38.gtf",
    "target_key": "GAPDH",
    "target_type": "gene",
    "primer_args": {
      "PRIMER_NUM_RETURN": 50,
      "PRIMER_OPT_TM": 60.0,
      "PRIMER_MIN_TM": 57.0,
      "PRIMER_MAX_TM": 63.0,
      "PRIMER_PRODUCT_SIZE_RANGE": "80-200"
    },
    "max_amplicon": 2000
  }
}
```

**Progress response** (sidecar → Electron, during long operations):
```json
{
  "id": 1,
  "status": "progress",
  "message": "Designing primers (junction mode)...",
  "pct": 30
}
```

**Final response** (sidecar → Electron):
```json
{
  "id": 1,
  "status": "ok",
  "result": {
    "filtered_pairs": [
      {
        "pair_name": "9746.1:45-199",
        "forward_seq": "TTCGTCGAAGGACTGCAGAC",
        "reverse_seq": "TGCAGTGCTTTCGAGACCAT",
        "forward_tm": 60.0,
        "reverse_tm": 60.0,
        "forward_gc": 55.0,
        "reverse_gc": 50.0,
        "product_size": 281
      }
    ],
    "warnings": [],
    "pdf_bytes_base64": "JVBERi0xLjQKJe...",
    "pdf_filename": "brimer-plast_GAPDH_2026-06-06.pdf"
  }
}
```

**Error response**:
```json
{
  "id": 1,
  "status": "error",
  "message": "No conserved exon chains found for gene 'GAPDH'"
}
```

Implementation sketch:
```python
import sys, json, base64, logging
from dataclasses import asdict
from pathlib import Path
from brimer_plast.pipeline import run_pipeline
from brimer_plast.pdf_report import build_pdf_report

def send(obj: dict) -> None:
    json.dump(obj, sys.stdout, default=str)
    sys.stdout.write("\n")
    sys.stdout.flush()

for line in sys.stdin:
    req = json.loads(line)
    rid = req["id"]
    
    if req["command"] == "run_pipeline":
        try:
            send({"id": rid, "status": "progress", "message": "Parsing annotations...", "pct": 10})
            result = run_pipeline(**req["params"])
            
            send({"id": rid, "status": "progress", "message": "Generating PDF report...", "pct": 90})
            pdf_bytes = build_pdf_report(result)  # needs to accept PipelineResult → bytes
            
            send({
                "id": rid,
                "status": "ok",
                "result": {
                    "filtered_pairs": [asdict(p) for p in result.filtered_pairs],
                    "warnings": result.warnings,
                    "pdf_bytes_base64": base64.b64encode(pdf_bytes).decode(),
                    "pdf_filename": f"brimer-plast_{req['params']['target_key']}_{datetime.now():%Y-%m-%d}.pdf",
                }
            })
        except Exception as e:
            send({"id": rid, "status": "error", "message": str(e)})
    else:
        send({"id": rid, "status": "error", "message": f"Unknown command: {req['command']}"})
```

**Decision points**:
- `build_pdf_report()` currently writes to a file path. Do we refactor it to
  return bytes? Or have the sidecar write to a temp path and return the path?
  → **Recommendation**: refactor to return `bytes` so the Electron app can
  trigger a browser download without touching the filesystem.

- Progress reporting: `run_pipeline()` currently has no progress hook. We
  would need to add one (a callback or generator). Is this worth the
  refactor, or do we accept a spinner with no percentage?
  → **Recommendation**: start with no progress (just a spinner), add it
  later if the 1–3 minute wait bothers users.

### 1.2 — Refactor `build_pdf_report` to return bytes

Current signature (from `pdf_report.py`):
```python
def build_pdf_report(
    result: PipelineResult,
    output_path: str | Path,
    ...
) -> None:
```

New signature:
```python
def build_pdf_report(
    result: PipelineResult,
    output_path: str | Path | None = None,
    ...
) -> bytes:
```
If `output_path` is None, render to a `BytesIO` buffer and return the bytes.
If `output_path` is given, write to file (backward compatible).

### 1.3 — PyInstaller bundle config

Create `pyinstaller.spec` or just use CLI flags. CI will run:

```bash
# On each platform (macOS, Windows):
pip install pyinstaller
pip install -e .                        # installs brimer-plast + deps
pyinstaller \
  --onefile \
  --name pybrimer \
  --add-binary tntblast;. \
  sidecar.py
```

This produces:
- macOS: `dist/pybrimer` (a Unix executable, actually a compressed bundle)
- Windows: `dist/pybrimer.exe`

The `tntblast` binary is placed next to the Python bundle, accessible via
`__file__` relative path from within the bundle.

---

## Phase 2 — Electron app

### 2.1 — `electron/package.json`

```json
{
  "name": "brimer-plast",
  "version": "0.1.0",
  "main": "main.js",
  "scripts": {
    "start": "electron .",
    "build:mac": "electron-builder --mac",
    "build:win": "electron-builder --win"
  },
  "devDependencies": {
    "electron": "^33.0.0",
    "electron-builder": "^25.0.0"
  },
  "build": {
    "appId": "com.brimer-plast",
    "productName": "Brimer-PLAST",
    "mac": {
      "target": "dmg",
      "icon": "build/icon.icns"
    },
    "win": {
      "target": "nsis",
      "icon": "build/icon.ico"
    },
    "extraResources": [
      {
        "from": "../dist/pybrimer",
        "to": "pybrimer"
      }
    ]
  }
}
```

### 2.2 — `electron/main.js` (~100 lines)

Responsibilities:
1. Create a `BrowserWindow` (900×700, no devtools in production)
2. On app ready, locate and spawn the Python sidecar:
   - Dev mode: `python sidecar.py` (look for `python3` or `python` on PATH)
   - Production: `process.resourcesPath + '/pybrimer'` (or `.exe` on Windows)
3. Pipe JSON:
   - Renderer sends IPC message → main process writes JSON to sidecar stdin
   - Sidecar stdout lines → main process parses JSON → sends IPC to renderer
4. On app quit, kill the sidecar process

**Key decision**: How does the renderer know when the sidecar is ready?
The sidecar has no startup handshake — it just blocks on stdin. First write
wakes it up. This is fine.

**Edge case**: If the sidecar crashes, `main.js` should detect the process
exit, log the error, and show a dialog. It should NOT auto-restart (that
could cause confusing duplicate runs).

### 2.3 — `electron/preload.js` (~20 lines)

```javascript
const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('api', {
  runPipeline: (params) => ipcRenderer.invoke('run-pipeline', params),
});
```

### 2.4 — `electron/renderer/index.html` (~350 lines)

Vanilla HTML + JS. No React, no build step.

HTML structure:
```
<form id="params">
  <!-- File inputs -->
  <label>Genome FASTA <input type="file" id="genome" accept=".fna,.fa,.fasta"></label>
  <label>GTF Annotation <input type="file" id="annotations" accept=".gtf"></label>

  <!-- Text inputs -->
  <label>Target gene(s) <textarea id="targets" rows="3" placeholder="GAPDH&#10;ACTB"></textarea></label>

  <!-- Collapsible advanced params -->
  <details>
    <summary>Primer parameters</summary>
    <label>Product min <input type="number" id="product-min" value="80"></label>
    <label>Product max <input type="number" id="product-max" value="200"></label>
    <!-- ... Tm, GC%, primer length ... -->
  </details>

  <button id="run-btn" type="submit">Run</button>
</form>

<div id="progress" style="display:none">
  <progress id="progress-bar"></progress>
  <span id="progress-text"></span>
</div>

<div id="results" style="display:none">
  <table id="pairs-table">
    <thead><tr>
      <th>Name</th>
      <th>Forward (5'→3')</th>
      <th>Reverse (5'→3')</th>
      <th>Tm (°C)</th>
      <th>GC%</th>
      <th>Size (bp)</th>
    </tr></thead>
    <tbody id="pairs-body"></tbody>
  </table>
  <button id="download-pdf">Download PDF report</button>
</div>

<div id="error" style="display:none; color:red"></div>
```

JS logic:
- On submit: disable button, show progress, collect form values
- Call `window.api.runPipeline(params)` — returns a Promise
- On progress IPC: update progress bar
- On result: populate table rows, enable PDF download
- On error: show error message

### 2.5 — Icons

Electron builder needs icons:
- macOS: `build/icon.icns` (1024×1024)
- Windows: `build/icon.ico` (256×256)

Options:
- Hire a designer for a logo
- Use a simple text-based icon ("BP" in a colored square)
- Use a public-domain DNA helix icon

---

## Phase 3 — CI/CD

### 3.1 — `.github/workflows/desktop.yml`

```yaml
name: Build desktop app

on:
  push:
    tags: ['v*']        # Only build on version tags
  workflow_dispatch:     # Or manual trigger

jobs:
  tntblast:
    strategy:
      matrix:
        include:
          - os: ubuntu-24.04   # Cross-compile for Windows via MinGW
            target: windows
            arch: x86_64
          - os: macos-13        # Intel Mac
            target: macos
            arch: x86_64
          - os: macos-14        # ARM Mac
            target: macos
            arch: arm64
    runs-on: ${{ matrix.os }}
    steps:
      - name: Download tnBLAST source
        run: |
          curl -sL -o tnblast.tar.gz \
            https://github.com/jgans/thermonucleotideBLAST/archive/refs/tags/v2.77.tar.gz
          mkdir tnblast-src && tar xzf tnblast.tar.gz -C tnblast-src --strip-components=1

      - name: Build (macOS)
        if: runner.os == 'macOS'
        run: |
          cd tnblast-src
          make CC=clang++ FLAGS="-O3 -Wall -std=c++14" LIBS="-lm -lz" -j$(sysctl -n hw.logicalcpu)

      - name: Setup MSYS2 (Windows cross-compile)
        if: runner.os == 'Linux'
        uses: msys2/setup-msys2@v2
        with:
          update: true
          install: >-
            mingw-w64-x86_64-gcc
            mingw-w64-x86_64-zlib

      - name: Build (Windows cross-compile)
        if: runner.os == 'Linux'
        shell: msys2 {0}
        run: |
          cd tnblast-src
          mingw32-make CC=g++ FLAGS="-O3 -Wall -std=c++14" LIBS="-lm -lz"

      - uses: actions/upload-artifact@v4
        with:
          name: tntblast-${{ matrix.target }}-${{ matrix.arch }}
          path: tnblast-src/tntblast*

  pybrimer:
    needs: tntblast
    strategy:
      matrix:
        include:
          - os: macos-13
            target: macos-x86_64
          - os: macos-14
            target: macos-arm64
          - os: windows-latest
            target: windows-x86_64
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: '3.12'

      - name: Install Python dependencies
        run: pip install . pyinstaller

      - name: Download tnBLAST binary
        uses: actions/download-artifact@v4
        with:
          name: tntblast-${{ matrix.target }}
          path: dist/bin/

      - name: Bundle with PyInstaller
        run: |
          pyinstaller \
            --onefile \
            --name pybrimer \
            --add-binary "dist/bin/tntblast*;." \
            sidecar.py

      - uses: actions/upload-artifact@v4
        with:
          name: pybrimer-${{ matrix.target }}
          path: dist/pybrimer*

  electron:
    needs: pybrimer
    strategy:
      matrix:
        include:
          - os: macos-13
            target: macos-x86_64
            builder-flag: --mac
          - os: macos-14
            target: macos-arm64
            builder-flag: --mac
          - os: windows-latest
            target: windows-x86_64
            builder-flag: --win
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node
        uses: actions/setup-node@v4
        with:
          node-version: 22

      - name: Download PyInstaller bundle
        uses: actions/download-artifact@v4
        with:
          name: pybrimer-${{ matrix.target }}
          path: dist/

      - name: Install npm dependencies
        run: npm install
        working-directory: electron

      - name: Build Electron installer
        run: npx electron-builder ${{ matrix.builder-flag }}
        working-directory: electron

      - uses: actions/upload-artifact@v4
        with:
          name: brimer-plast-${{ matrix.target }}
          path: electron/dist/*.dmg
          if-no-files-found: error
```

### 3.2 — Code signing

Without signing:

| Platform | User sees | Workaround |
|---|---|---|
| macOS | "Brimer-PLAST.app is damaged and can't be opened" | `xattr -d com.apple.quarantine /Applications/Brimer-PLAST.app` |
| Windows | "Windows protected your PC" (SmartScreen) | Click "Run anyway" |

For early releases, you can distribute without signing and tell users how to
bypass the warnings. For wider distribution, budget for:
- Apple Developer account: $99/year
- Windows code signing certificate: ~$200/year (or use Azure Key Vault + eSigner)

### 3.3 — CI testing

Add a post-build smoke test job that:
1. Downloads the Electron installer
2. Installs it on the runner
3. Runs the app with `--no-sandbox` and a small test genome
4. Verifies the PDF output is non-empty

Tests are on the C. elegans genome already in `tests/fixtures/ce11/` (~42 MB).

---

## Decision log

| # | Question | Decision | Rationale |
|---|---|---|---|
| 1 | Where does tnBLAST compile? | CI downloads upstream source, builds, discards | No submodule, no fork, no vendored code to maintain |
| 2 | tnBLAST OpenMP? | Disabled | Single-threaded is fast enough; avoids macOS libomp dependency |
| 3 | Progress reporting from run_pipeline? | Not yet — just a spinner | Adding a progress callback to run_pipeline() is a separate refactor |
| 4 | PDF generation: file vs bytes? | Return bytes | Cleaner for Electron download flow |
| 5 | Electron with vanilla HTML or React? | Vanilla HTML | No build step, no npm deps beyond electron itself |
| 6 | GitHub Actions or self-hosted? | GitHub Actions free tier | No infrastructure to maintain |
| 7 | Code signing for first release? | Skip | Distribute unsigned with instructions; add signing later |
| 8 | macOS ARM or Intel only? | Both | GitHub provides ARM runners for free; no reason to exclude ARM users |

---

## Effort estimate (updated)

| Phase | Days | Risk factor |
|---|---|---|
| 0.1 — tnBLAST on macOS | 0.5 | Low — standard C++14, no platform code |
| 0.2 — tnBLAST on Windows (MinGW) | 1 | Medium — MinGW setup, `getopt_long` availability |
| 0.3 — primer3-py on ARM | 0.5 | Medium — libprimer3 Makefile with Clang on ARM |
| 1.1 — sidecar.py | 0.5 | Trivial |
| 1.2 — Refactor build_pdf_report | 0.5 | Low |
| 1.3 — PyInstaller bundle | 0.5 | Low — standard tooling |
| 2.1–2.4 — Electron app | 2 | Medium — HTML GUI design, IPC debugging |
| 3.1 — CI workflow | 1.5 | Medium — GitHub Actions matrix debugging |
| 3.2 — Testing on clean OS | 1 | Low — just run the installer |
| **Total** | **~8 days** | |

The biggest risk is **0.2 (Windows MinGW + getopt_long)** — if MinGW's
implementation of `getopt_long` behaves differently, we need a fallback.
Everything else is low-risk.

---

## Next steps

1. **Start with Phase 0.1**: Create a one-shot CI job that downloads tnBLAST
   source and compiles it on `macos-14`. See if it works. Takes 15 minutes.
2. **Then 0.2**: Same for `windows-latest` with MSYS2.
3. **Then 0.3**: Add `pip install primer3-py` to the macOS ARM job and check
   the exit code.
4. **Then phases 1–3 in order**.

None of Phase 0 blocks anything in Phase 1 or 2 — you could write sidecar.py
and the Electron app on your current machine right now, without any cross-
compiled binaries. The sidecar just needs `python sidecar.py` for development.