#!/bin/bash
# ================================
# run_analysis.sh
# Usage: ./run_analysis.sh <GFF_FILE> [FASTA_FILE] [OUT_DIR] [--strict]
# ================================

GFF_FILE="$1"
FASTA_FILE=""
OUT_DIR=""
STRICT_FLAG=""

# Determine strict flag (if last argument is --strict)
if [ "${@: -1}" == "--strict" ]; then
    STRICT_FLAG="--strict"
fi

# Determine FASTA and OUT_DIR based on arguments
if [ -z "$3" ]; then
    # Only 2 arguments → GFF and OUT_DIR
    OUT_DIR="$2"
else
    # 3 arguments → FASTA provided
    FASTA_FILE="$2"
    OUT_DIR="$3"
fi


# Create output folder if missing
mkdir -p "$OUT_DIR"

# Get directory from GFF for Docker mount
INPUT_DIR="$(cd "$(dirname "$GFF_FILE")" && pwd)"

# Print info
echo "GFF file: $GFF_FILE"
echo "FASTA file: $FASTA_FILE"
echo "Output directory: $OUT_DIR"
echo "Mounted input directory: $INPUT_DIR"
echo "Passing strict flag: $STRICT_FLAG"

# Build Docker command
DOCKER_CMD=(docker run --rm -v "$INPUT_DIR:/work" -v "$OUT_DIR:/output" -w /work gene-summariser -g "$(basename "$GFF_FILE")" --outdir "/output" "$STRICT_FLAG")

# Add FASTA if provided
if [ -n "$FASTA_FILE" ]; then
    DOCKER_CMD+=(-f "$(basename "$FASTA_FILE")")
fi

# Add strict flag if provided
if [ -n "$STRICT_FLAG" ]; then
    DOCKER_CMD+=("$STRICT_FLAG")
fi

# Run Docker
"${DOCKER_CMD[@]}"