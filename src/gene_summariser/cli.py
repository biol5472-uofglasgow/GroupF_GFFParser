import argparse
import os
import sys
from pathlib import Path

from gene_summariser.figures import (
    CDSLengthHistogram,
    ExonCountHistogram,
    FlaggedBarChart,
    PieChart,
)
from gene_summariser.metrics import MetricsCalculator

# importing required modules from the gene_summariser package
from gene_summariser.parser import ParserGFF
from gene_summariser.qc import QCChecker
from gene_summariser.writer import OutputWriter


def main() -> None:
    """
    Main entry point for the CLI
    This function :
    - parses command line arguments
    - validates the user input
    - Orchestrates parsing, QC checks, and output writing
    - Handles errors and strict-mode behaviour if required
    """

    # Argument parsing
    parser = argparse.ArgumentParser(
        prog="QC CHECK ON GFF",
        description="""This program takes a gff file and performs necessary QC checks on
        it to ensure that everything is correct.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required input GFF file
    parser.add_argument(
        "-g", "--gff", required=True, help="please write the path of the gffFile here."
    )

    # Optional genome FASTA file
    parser.add_argument(
        "-f",
        "--fasta",
        help="Optional ,please write the path of the genome fasta file here.",
    )

    # Output options
    parser.add_argument(
        "-o",
        "--output",
        default="qc_report.txt",
        help="Output QC report file (default: qc_report.txt)",
    )

    parser.add_argument(
        "--log", default="qc.log", help="Log file for detailed execution info"
    )

    # fail if any strict mode is detected
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail execution if any QC warning is detected",
    )

    # Output format option
    parser.add_argument(
        "--format",
        choices=["text", "csv", "json"],
        default="text",
        help="Output format for QC report",
    )
    parser.add_argument(
        "--outdir", default="results", help="Output directory for all generated files"
    )

    args = parser.parse_args()

    # checking file
    if not os.path.isfile(args.gff):
        print(f"GFF file not found: {args.gff}")
        SystemExit(1)

    if args.fasta and not os.path.isfile(args.fasta):
        print(f"FASTA file not found: {args.fasta}")
        SystemExit(1)

    # checking file type
    if not args.gff.endswith((".gff", ".gff3")):
        print("Input file is not a GFF/GFF3 file")
        SystemExit(1)

    try:
        print("Parsing GFF file...")
        parser_gff = ParserGFF(args.gff)
        transcripts = parser_gff.parse_transcripts()
        print(f" Parsed{len(transcripts)} transcripts")

        print("Running QC checks...")
        qc_checker = QCChecker(fasta_file=args.fasta)

        print("Calculating transcripts metrics and qc flags")
        calculator = MetricsCalculator(qc_checker)
        summaries = calculator.calculate_summaries(transcripts)

        print("Writing output")
        outdir = Path(args.outdir)
        writer = OutputWriter(outdir)
        summary_path = writer.write_transcript_summary(summaries)
        writer.write_provenance(
            input_file=Path(args.gff),
            parameters={
                "fasta": args.fasta,
                "strict": args.strict,
            },
        )
        writer.write_qc_flags_gff3(transcripts, summaries)
        writer.write_qc_flags_bed(transcripts, summaries)

        pie_chart = PieChart(summaries, output_dir=outdir)
        pie_chart.generate_pie_chart()

        exon_histogram = ExonCountHistogram(summaries, output_dir=outdir)
        exon_histogram.generate_histogram()

        flagged_bar_chart = FlaggedBarChart(summaries, output_dir=outdir)
        flagged_bar_chart.generate_bar_plot()

        cds_length_histogram = CDSLengthHistogram(summaries, output_dir=outdir)
        cds_length_histogram.generate_histogram()

        print(f"  Transcript summary written to: {summary_path}")

        if args.strict:
            n_flagged = sum(1 for s in summaries if s.flags)
            if n_flagged > 0:
                print(
                    f"ERROR: mode failed — {n_flagged} transcripts have QC issues.\n"
                    "Inspect 'transcript_summary.tsv' or 'qc_flags.gff3' for details.\n"
                    "Run without --strict to generate reports without failing.",
                    file=sys.stderr,
                )
                sys.exit(2)

        print("QC checks completed")

    except Exception as e:
        print(f"Error:{e}")
        SystemExit(1)


if __name__ == "__main__":
    main()
