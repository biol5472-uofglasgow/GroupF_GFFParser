"""Generate HTML report from gene-summariser results."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------------
# Figure generation (mirrors the logic in figures.py, but produces files
# relative to the report so the HTML can reference them with simple paths)
# ---------------------------------------------------------------------------

def _ensure_figures_dir(figures_dir: Path) -> Path:
    figures_dir.mkdir(parents=True, exist_ok=True)
    return figures_dir


def _generate_pie_chart(flag_counts: dict[str, int], figures_dir: Path) -> str:
    """Pie chart of flag distribution. Returns the relative src path."""
    filename = "pie_chart.png"
    labels = list(flag_counts.keys())
    values = list(flag_counts.values())

    plt.figure(figsize=(8, 8))
    plt.pie(values, startangle=90)
    plt.legend(labels, bbox_to_anchor=(1, 0.5))
    plt.title("Flag Distribution")
    plt.axis("equal")
    plt.savefig(figures_dir / filename, bbox_inches="tight")
    plt.close()
    return f"figures/{filename}"


def _generate_flagged_bar(flagged: int, unflagged: int, figures_dir: Path) -> str:
    """Bar chart of flagged vs unflagged. Returns the relative src path."""
    filename = "flagged_vs_unflagged_bar.png"
    group  = ["flagged", "unflagged"]
    counts = [flagged, unflagged]
    labels = ["Flagged", "Unflagged"]
    colors = ["#409ab6", "#ce6868"]

    plt.figure(figsize=(8, 8))
    plt.bar(group, counts, label=labels, color=colors)
    plt.xlabel("Transcript Status")
    plt.ylabel("Number of transcripts")
    plt.title("Flagged vs Unflagged Transcripts")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / filename)
    plt.close()
    return f"figures/{filename}"


def _generate_exon_histogram(exon_counts: list[int], figures_dir: Path) -> str:
    """Histogram of exon counts. Returns the relative src path."""
    filename = "exon_count_distribution.png"

    plt.figure(figsize=(8, 8))
    plt.hist(
        exon_counts,
        bins=range(1, max(exon_counts) + 2),
        edgecolor="black",
        alpha=0.7,
    )
    plt.xlabel("Number of Exons")
    plt.ylabel("Number of Transcripts")
    plt.title("Exon Count Distribution")
    plt.xticks(range(1, max(exon_counts) + 1))
    plt.grid(axis="y", alpha=0.75)
    plt.savefig(figures_dir / filename)
    plt.close()
    return f"figures/{filename}"


def _generate_cds_histogram(cds_lengths: list[int | float], figures_dir: Path) -> str | None:
    """Histogram of CDS lengths. Returns the relative src path, or None if no data."""
    if not cds_lengths:
        return None

    filename = "cds_length_distribution.png"

    plt.figure(figsize=(8, 8))
    plt.hist(cds_lengths, edgecolor="black", alpha=0.7)
    plt.xlabel("CDS Length")
    plt.ylabel("Number of Transcripts")
    plt.title("CDS Length Distribution")
    plt.grid(axis="y", alpha=0.75)
    plt.savefig(figures_dir / filename)
    plt.close()
    return f"figures/{filename}"


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def generate_html_report(results_dir: Path, output_file: Path) -> None:
    """Generate HTML report with visualizations."""

    # Load data
    summary_df = pd.read_csv(results_dir / "transcript_summary.tsv", sep="\t")

    with open(results_dir / "run.json") as f:
        provenance = json.load(f)

    # --------------- statistics ------------------------------------------------
    total_transcripts = len(summary_df)
    flagged_count     = len(summary_df[summary_df["flags"] != ""])
    unflagged_count   = total_transcripts - flagged_count
    with_cds          = summary_df["has_cds"].sum()

    # --------------- flag counts -----------------------------------------------
    flag_counts: dict[str, int] = {"no_flag": 0}
    for flags in summary_df["flags"]:
        if flags:
            for flag in flags.split(","):
                flag_counts[flag] = flag_counts.get(flag, 0) + 1
        else:
            flag_counts["no_flag"] += 1

    # --------------- exon / CDS lists ------------------------------------------
    exon_counts  = summary_df["n_exons"].tolist()
    cds_lengths  = [
        v for v in summary_df["total_cds_length"]
        if pd.notna(v) and v > 0
    ]

    # --------------- generate figures ------------------------------------------
    # Figures are saved next to the output HTML so relative paths work
    figures_dir = _ensure_figures_dir(output_file.parent / "figures")

    pie_src     = _generate_pie_chart(flag_counts, figures_dir)
    bar_src     = _generate_flagged_bar(flagged_count, unflagged_count, figures_dir)
    exon_src    = _generate_exon_histogram(exon_counts, figures_dir)
    cds_src     = _generate_cds_histogram(cds_lengths, figures_dir)

    # --------------- build figure HTML snippet ---------------------------------
    def _img(src: str | None, alt: str) -> str:
        if src is None:
            return ""
        return (
            f'<div class="figure">'
            f'<img src="{src}" alt="{alt}">'
            f'<p class="figure-caption">{alt}</p>'
            f'</div>'
        )

    figures_html = (
        _img(pie_src,  "Flag Distribution") +
        _img(bar_src,  "Flagged vs Unflagged Transcripts") +
        _img(exon_src, "Exon Count Distribution") +
        _img(cds_src,  "CDS Length Distribution")
    )

    # --------------- flag summary rows -----------------------------------------
    flag_rows = ""
    for flag, count in sorted(flag_counts.items(), key=lambda x: x[1], reverse=True):
        pct = count / total_transcripts * 100
        flag_rows += f"        <tr><td>{flag}</td><td>{count}</td><td>{pct:.1f}%</td></tr>\n"

    # --------------- flagged transcript rows -----------------------------------
    transcript_rows = ""
    for _, row in summary_df[summary_df["flags"] != ""].head(50).iterrows():
        flags_html = " ".join(
            f'<span class="flag">{f}</span>' for f in row["flags"].split(",") if f
        )
        transcript_rows += (
            f"        <tr>"
            f"<td>{row['gene_id']}</td>"
            f"<td>{row['transcript_id']}</td>"
            f"<td>{row['n_exons']}</td>"
            f"<td>{'Yes' if row['has_cds'] else 'No'}</td>"
            f"<td>{flags_html}</td>"
            f"</tr>\n"
        )

    # --------------- assemble HTML ---------------------------------------------
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

        .figures {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 30px; margin: 20px 0; }}
        .figure {{ text-align: center; }}
        .figure img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }}
        .figure-caption {{ font-size: 0.9em; color: #555; margin-top: 6px; }}
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
            <div class="stat-value">{flagged_count}</div>
            <div>Flagged Transcripts</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{with_cds}</div>
            <div>With CDS</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{flagged_count/total_transcripts*100:.1f}%</div>
            <div>Flag Rate</div>
        </div>
    </div>

    <h2>QC Flag Summary</h2>
    <table>
        <tr><th>Flag</th><th>Count</th><th>Percentage</th></tr>
{flag_rows}    </table>

    <h2>Visualizations</h2>
    <div class="figures">
{figures_html}
    </div>

    <h2>Flagged Transcripts (First 50)</h2>
    <table>
        <tr><th>Gene</th><th>Transcript</th><th>Exons</th><th>CDS</th><th>Flags</th></tr>
{transcript_rows}    </table>

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