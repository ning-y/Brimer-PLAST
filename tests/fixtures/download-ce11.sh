#!/usr/bin/env bash
# Download C. elegans (WBcel235) genome and annotations for integration tests.
#
# Downloads to tests/fixtures/ce11/  (about 42 MB total).
# Run from the repository root.

set -euo pipefail

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235"
OUT_DIR="tests/fixtures/ce11"

mkdir -p "$OUT_DIR"

echo "Downloading C. elegans genome (WBcel235)..."
curl -sL "$BASE_URL/GCF_000002985.6_WBcel235_genomic.fna.gz" \
    -o "$OUT_DIR/genome.fna.gz" &
curl -sL "$BASE_URL/GCF_000002985.6_WBcel235_genomic.gtf.gz" \
    -o "$OUT_DIR/annotations.gtf.gz" &

wait

echo "Decompressing..."
gunzip -f "$OUT_DIR/genome.fna.gz"
gunzip -f "$OUT_DIR/annotations.gtf.gz"

echo "Done."
ls -lh "$OUT_DIR/"