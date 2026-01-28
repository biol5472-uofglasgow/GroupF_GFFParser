#!/bin/bash
# Simple wrapper script for gene-summariser

set -e  # Exit on error

# Default values
GFF_FILE=""
FASTA_FILE=""
OUTPUT_DIR="results"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --gff)
      GFF_FILE="$2"
      shift 2
      ;;
    --fasta)
      FASTA_FILE="$2"
      shift 2
      ;;
    --outdir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Validate
if [ -z "$GFF_FILE" ]; then
  echo "Error: --gff required"
  exit 1
fi

# Run
echo "Running gene-summariser..."
if [ -n "$FASTA_FILE" ]; then
  gene-summariser --gff "$GFF_FILE" --fasta "$FASTA_FILE" --outdir "$OUTPUT_DIR"
else
  gene-summariser --gff "$GFF_FILE" --outdir "$OUTPUT_DIR"
fi

echo "Done! Results in $OUTPUT_DIR"