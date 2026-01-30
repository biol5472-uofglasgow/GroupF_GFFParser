#!/bin/bash
# Usage: ./docker-run-test.sh <gff_file> <fasta_file> <outdir>
# Must be run from the project root

GFF_FILE="$1"
FASTA_FILE="$2"
OUT_DIR="$3"

# Set default output folder
if [ -z "$OUT_DIR" ]; then
  OUT_DIR="results"
fi

mkdir -p "$OUT_DIR"

docker run --rm \
  -v "$PWD":/work \
  -w /work \
  gene-summariser:latest \
  --gff "$GFF_FILE" \
  --fasta "$FASTA_FILE" \
  --outdir "$OUT_DIR"
