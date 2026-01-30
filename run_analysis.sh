#!/bin/bash
# ================================
# run_analysis.sh
# Usage: ./run_analysis.sh <GFF_FILE> <FASTA_FILE> [OUT_DIR] [--strict]
# ================================

GFF_FILE="$1"
FASTA_FILE="$2"
OUT_DIR="${3:-results}"  # default to 'results' if not provided

# Check if last argument is --strict
STRICT_FLAG=""
if [ "${@: -1}" == "--strict" ]; then
    STRICT_FLAG="--strict"
fi

# Check GFF file exists
if [ ! -f "$GFF_FILE" ]; then
    echo "ERROR: GFF file not found: $GFF_FILE"
    exit 1
fi

# Check FASTA file exists, warn if not
if [ ! -f "$FASTA_FILE" ]; then
    echo "WARNING: FASTA file not found: $FASTA_FILE"
fi

# Make sure output directory exists
mkdir -p "$OUT_DIR"

# Get directory of GFF file (for Docker mount)
INPUT_DIR="$(dirname "$GFF_FILE")"

# Print info
echo "GFF file: $GFF_FILE"
echo "FASTA file: $FASTA_FILE"
echo "Output directory: $OUT_DIR"
echo "Mounted input directory: $INPUT_DIR"
echo "Passing strict flag: $STRICT_FLAG"

# Run Docker
docker run --rm -v "$INPUT_DIR:/work" -w /work gene-summariser \
    -g "$(basename "$GFF_FILE")" \
    -f "$(basename "$FASTA_FILE")" \
    --outdir "$OUT_DIR" $STRICT_FLAG