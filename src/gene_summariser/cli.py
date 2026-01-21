import argparse
import os 
import sys

def main():
    parser = argparse.ArgumentParser(prog='QC CHECK ON GFF',description='This program takes a gff file and performs necessary QC checks on it to ensure that everything is correct.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-g","--gff", required=True, help='please write the path of the gffFile here.')
    parser.add_argument("-f","--fasta",help='Optional ,please write the path of the genome fasta file here.')
    parser.add_argument("-o", "--output",default="qc_report.txt",help="Output QC report file (default: qc_report.txt)")
    parser.add_argument("--log",default="qc.log",help="Log file for detailed execution info")
    parser.add_argument("--strict",action="store_true",help="Fail execution if any QC warning is detected")
    parser.add_argument("--format",choices=["text", "csv", "json"],default="text",help="Output format for QC report")

    return parser.parse_args()

# def validate_inputs(args):
#     if not os.path.isfile(args.gff):
#         sys.exit(f"GFF file not found: {args.gff}")

#     if args.fasta and not os.path.isfile(args.fasta):
#         sys.exit(f"FASTA file not found: {args.fasta}")

#     if not args.gff.endswith((".gff", ".gff3")):
#         sys.exit("Input file is not a GFF/GFF3 file")

    
