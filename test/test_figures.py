import os
from pathlib import Path

import pytest

from gene_summariser.figures import (
    CDSLengthHistogram,
    ExonCountHistogram,
    FlaggedBarChart,
    PieChart,
)
from gene_summariser.models import TranscriptSummary


@pytest.fixture
def transcripts():
    return [
        TranscriptSummary(
            transcript_id="tx1",
            gene_id="gene1",
            flags=["flag1", "flag2"],
            n_exons=2,
            has_cds=True,
        ),
        TranscriptSummary(
            transcript_id="tx2",
            gene_id="gene2",
            flags=["flag1"],
            n_exons=2,
            has_cds=True,
        ),
        TranscriptSummary(
            transcript_id="tx3",
            gene_id="gene3",
            flags=[],
            n_exons=1,
            has_cds=True,
        ),
    ]


def test_pie_chart_generation(
    transcripts: list[TranscriptSummary], tmp_path: Path
) -> None:
    """
    Test that pie chart generation works without error
    """
    output_file = tmp_path / "figures" / "pie_chart.png"
    figures_dir = tmp_path / "figures"
    piechart = PieChart(transcripts, output_dir=tmp_path)
    piechart.generate_pie_chart()

    assert figures_dir.exists()
    assert os.path.exists(output_file)


def test_exon_histogram_generation(
    transcripts: list[TranscriptSummary], tmp_path: Path
) -> None:
    """
    Test that exon histogram generation works without error
    """
    output_file = tmp_path / "figures" / "exon_count_distribution.png"
    figures_dir = tmp_path / "figures"
    exon_histogram = ExonCountHistogram(transcripts, output_dir=tmp_path)
    exon_histogram.generate_histogram()

    assert figures_dir.exists()
    assert os.path.exists(output_file)


def test_flagged_bar_chart_generation(
    transcripts: list[TranscriptSummary], tmp_path: Path
) -> None:
    """
    Test that flagged bar chart generation works without error
    """
    output_file = tmp_path / "figures" / "flagged_vs_unflagged_bar.png"
    figures_dir = tmp_path / "figures"
    flagged_bar_chart = FlaggedBarChart(transcripts, output_dir=tmp_path)
    flagged_bar_chart.generate_bar_plot()

    assert figures_dir.exists()
    assert os.path.exists(output_file)


def test_cds_length_histogram_generation(
    transcripts: list[TranscriptSummary], tmp_path: Path
) -> None:
    """
    Test that CDS length histogram generation works without error
    """
    output_file = tmp_path / "figures" / "cds_length_distribution.png"
    figures_dir = tmp_path / "figures"
    cds_length_histogram = CDSLengthHistogram(transcripts, output_dir=tmp_path)
    cds_length_histogram.generate_histogram()

    assert figures_dir.exists()
    assert os.path.exists(output_file)
