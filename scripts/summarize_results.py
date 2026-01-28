"""Command-line tool to summarize QC results."""

import argparse
import json
from pathlib import Path
from collections import Counter
import pandas as pd


def print_summary(results_dir: Path) -> None:
    """Print summary of QC results."""
    
    # Load data
    df = pd.read_csv(results_dir / "transcript_summary.tsv", sep="\t")
    
    with open(results_dir / "run.json") as f:
        prov = json.load(f)
    
    # Print header
    print("=" * 70)
    print(f"Gene Model QC Summary Report")
    print("=" * 70)
    print()
    
    # Provenance
    print(f"Input File: {prov['inputs']['gff3_file']}")
    print(f"Analysis Date: {prov['timestamp']}")
    print(f"Tool Version: {prov['version']}")
    print()
    
    # Basic stats
    total = len(df)
    flagged = len(df[df['flags'] != ''])
    with_cds = df['has_cds'].sum()
    
    print("Overview:")
    print(f"  Total Transcripts:   {total:>10,}")
    print(f"  Flagged Transcripts: {flagged:>10,} ({flagged/total*100:.1f}%)")
    print(f"  With CDS:            {with_cds:>10,} ({with_cds/total*100:.1f}%)")
    print(f"  Without CDS:         {total-with_cds:>10,} ({(total-with_cds)/total*100:.1f}%)")
    print()
    
    # Exon stats
    print("Exon Statistics:")
    print(f"  Mean exons per transcript: {df['n_exons'].mean():.1f}")
    print(f"  Median exons:              {df['n_exons'].median():.0f}")
    print(f"  Max exons:                 {df['n_exons'].max()}")
    print(f"  Single-exon transcripts:   {(df['n_exons'] == 1).sum()}")
    print()
    
    # Flag breakdown
    flag_counts = Counter()
    for flags in df['flags']:
        if flags:
            for flag in flags.split(","):
                flag_counts[flag] += 1
    
    if flag_counts:
        print("QC Flags:")
        for flag, count in flag_counts.most_common():
            print(f"  {flag:<30} {count:>6} ({count/total*100:>5.1f}%)")
    else:
        print("QC Flags: None (all transcripts passed QC!)")
    print()
    
    # Critical issues
    critical_flags = ['no_exons', 'overlapping_exons', 'cds_phase_inconsistent']
    critical_count = sum(flag_counts[f] for f in critical_flags if f in flag_counts)
    
    if critical_count > 0:
        print(f"⚠️  WARNING: {critical_count} transcripts with CRITICAL issues!")
        for flag in critical_flags:
            if flag in flag_counts:
                print(f"    - {flag}: {flag_counts[flag]}")
    else:
        print("✓ No critical issues found")
    
    print()
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize gene-summariser results"
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="Directory containing results"
    )
    args = parser.parse_args()
    
    if not (args.results_dir / "transcript_summary.tsv").exists():
        print(f"Error: No results found in {args.results_dir}")
        print("Run gene-summariser first!")
        return 1
    
    print_summary(args.results_dir)
    return 0


if __name__ == "__main__":
    exit(main())