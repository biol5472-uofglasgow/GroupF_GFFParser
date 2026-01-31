
import os
import sys


def validate_inputs(args):
    if not os.path.isfile(args.gff):
        sys.exit(f"GFF file not found: {args.gff}")

    if args.fasta and not os.path.isfile(args.fasta):
        sys.exit(f"FASTA file not found: {args.fasta}")

    if not args.gff.endswith((".gff", ".gff3")):
        sys.exit("Input file is not a GFF/GFF3 file")
