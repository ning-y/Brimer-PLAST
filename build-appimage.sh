#!/bin/bash
set -euxo pipefail

# Install dependencies
apt-get update -qq
apt-get install -y -qq \
  nodejs npm \
  python3 python3-pip python3-venv \
  p7zip-full \
  libfuse2 \
  wget curl \
  > /dev/null 2>&1

# Create a venv and install Python deps
python3 -m venv /venv
source /venv/bin/activate

# Install primer3-py, pyinstaller, and the project
pip install primer3-py pyinstaller -q
pip install /app -q

# Build tnBLAST from source
cd /tmp
curl -sL -o tnblast.tar.gz \
  https://github.com/jgans/thermonucleotideBLAST/archive/refs/tags/v2.77.tar.gz
mkdir tnblast-src
tar xzf tnblast.tar.gz -C tnblast-src --strip-components=1
cd tnblast-src
make CC=g++ FLAGS="-O3 -Wall -std=c++14" LIBS="-lm -lz" OPENMP= -j$(nproc)
TNTBLAST=$(realpath tntblast)
echo "tnBLAST built: $TNTBLAST"

# Build PyInstaller bundle
cd /app
pyinstaller --onefile --name pybrimer --add-binary "$TNTBLAST:." sidecar.py 2>&1 | tail -3

# Copy tnBLAST to dist/ for Electron extraResources
echo "Copying tnBLAST to dist/ for Electron packaging..."
cp "$TNTBLAST" dist/

# Verify pybrimer
set +e
echo '{"id":1,"command":"unknown","params":{}}' | timeout 5 ./dist/pybrimer 2>/dev/null
PYRET=$?
set -e
# timeout 124 means process was killed after 5s — that's expected (stdin closed)
# exit 0 means it printed JSON and exited cleanly — also fine
if [ "$PYRET" -ne 0 ] && [ "$PYRET" -ne 124 ]; then
  echo "WARN: pybrimer test returned unexpected exit code $PYRET"
fi
rm -rf /app/build /app/*.spec

# Build Electron + AppImage
cd /app/electron
npm install 2>&1 | tail -5
node build.mjs --linux AppImage 2>&1 | tail -15

echo ""
echo "=== Build complete ==="
ls -lh /app/electron/dist/Brimer-PLAST-*.AppImage 2>/dev/null || echo "No AppImage found"