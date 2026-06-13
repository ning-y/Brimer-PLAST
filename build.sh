#!/usr/bin/env bash
set -euo pipefail

# Build the Brimer-PLAST AppImage locally.
# Prerequisites: `nix develop` shell (provides Python, tnBLAST, nodejs).

echo "=== 1. Install Python build deps ==="
pip install pyinstaller -q
pip install -e . -q

echo "=== 2. Locate tnBLAST binary ==="
TNT=$(which tntblast)
echo "  tnBLAST: $TNT"

echo "=== 3. Build PyInstaller bundle ==="
pyinstaller \
  --onefile \
  --name pybrimer \
  --collect-all primer3 \
  --add-binary "$TNT:." \
  sidecar.py

echo "=== 3.5. Copy tnBLAST to dist/ for Electron resources ==="
cp "$TNT" dist/

echo "=== 4. Install npm deps ==="
cd electron
npm install --silent

echo "=== 5. Build AppImage ==="
npx electron-builder --linux AppImage

echo "=== Done ==="
ls -lh dist/*.AppImage 2>/dev/null || echo "Check electron/dist/ for output"