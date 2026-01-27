"""Generate HTML report from gene-summariser results."""

import argparse
import json
from pathlib import Path
import pandas as pd


def generate_html_report(results_dir: Path, output_file: Path) -> None:
    """Generate HTML report with visualizations."""
    
    # Load data
    summary_df = pd.read_csv(results_dir / "transcript_summary.tsv", sep="\t")
    
    with open(results_dir / "run.json") as f:
        provenance = json.load(f)
    
    # Calculate statistics
    total_transcripts = len(summary_df)
    flagged = len(summary_df[summary_df["flags"] != ""])
    with_cds = summary_df["has_cds"].sum()
    
    # Count flags
    flag_counts = {}
    for flags in summary_df["flags"]:
        if flags:
            for flag in flags.split(","):
                flag_counts[flag] = flag_counts.get(flag, 0) + 1
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Gene Model QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #4CAF50; color: white; padding: 20px; }}
        .stats {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 20px; margin: 20px 0; }}
        .stat-box {{ background: #f0f0f0; padding: 20px; text-align: center; }}
        .stat-value {{ font-size: 2em; font-weight: bold; color: #4CAF50; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #4CAF50; color: white; }}
        .flag {{ background: #ffebee; padding: 2px 6px; border-radius: 3px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Gene Model QC Report</h1>
        <p>Generated: {provenance['timestamp']}</p>
    </div>
    
    <div class="stats">
        <div class="stat-box">
            <div class="stat-value">{total_transcripts}</div>
            <div>Total Transcripts</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{flagged}</div>
            <div>Flagged Transcripts</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{with_cds}</div>
            <div>With CDS</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{flagged/total_transcripts*100:.1f}%</div>
            <div>Flag Rate</div>
        </div>
    </div>
    
    <h2>QC Flag Summary</h2>
    <table>
        <tr><th>Flag</th><th>Count</th><th>Percentage</th></tr>
"""
    
    for flag, count in sorted(flag_counts.items(), key=lambda x: x[1], reverse=True):
        pct = count / total_transcripts * 100
        html += f"<tr><td>{flag}</td><td>{count}</td><td>{pct:.1f}%</td></tr>"
    
    html += """
    </table>
    
    <h2>Flagged Transcripts (First 50)</h2>
    <table>
        <tr><th>Gene</th><th>Transcript</th><th>Exons</th><th>CDS</th><th>Flags</th></tr>
"""
    
    flagged_df = summary_df[summary_df["flags"] != ""].head(50)
    for _, row in flagged_df.iterrows():
        flags_html = " ".join([f'<span class="flag">{f}</span>' 
                              for f in row["flags"].split(",") if f])
        html += f"""
        <tr>
            <td>{row['gene_id']}</td>
            <td>{row['transcript_id']}</td>
            <td>{row['n_exons']}</td>
            <td>{'Yes' if row['has_cds'] else 'No'}</td>
            <td>{flags_html}</td>
        </tr>
"""
    
    html += f"""
    </table>
    
    <h2>Provenance</h2>
    <table>
        <tr><th>Tool</th><td>{provenance['tool']} v{provenance['version']}</td></tr>
        <tr><th>Input File</th><td>{provenance['inputs']['gff3_file']}</td></tr>
        <tr><th>Timestamp</th><td>{provenance['timestamp']}</td></tr>
        <tr><th>Parameters</th><td>{json.dumps(provenance['parameters'])}</td></tr>
    </table>
</body>
</html>
"""
    
    output_file.write_text(html)
    print(f"Report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate HTML QC report")
    parser.add_argument("results_dir", type=Path, help="Results directory")
    parser.add_argument("-o", "--output", type=Path, default=Path("qc_report.html"))
    args = parser.parse_args()
    
    generate_html_report(args.results_dir, args.output)
    print(f"\nOpen in browser: file://{args.output.absolute()}")


if __name__ == "__main__":
    main()