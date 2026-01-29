from pathlib import Path

import pytest

from gene_summariser.generate_report import generate_html_report
from gene_summariser.metrics import MetricsCalculator
from gene_summariser.parser import ParserGFF
from gene_summariser.qc import QCChecker
from gene_summariser.writer import OutputWriter


def project_pipeline(gff_file: str, fasta_file: str, outdir: Path) -> None:
    parser_gff = ParserGFF(gff_file)
    transcripts = parser_gff.parse_transcripts()

    qc_checker = QCChecker(fasta_file)
    calculator = MetricsCalculator(qc_checker)
    summaries = calculator.calculate_summaries(transcripts)

    writer = OutputWriter(outdir)
    summary_path = writer.write_transcript_summary(summaries)
    writer.write_provenance(
        input_file=Path(gff_file),
        parameters={"fasta": fasta_file},
    )
    writer.write_qc_flags_gff3(transcripts, summaries)
    writer.write_qc_flags_bed(transcripts, summaries)

    generate_html_report(outdir)


def test_integration(tmp_path: Path) -> None:
    gff_file = "test/fixtures/models.gff3"
    fasta_file = "test/fixtures/testfasta.fasta"
    outdir = tmp_path / "output"
    outdir.mkdir()

    project_pipeline(gff_file, fasta_file, outdir)

    summary_path = outdir / "transcript_summary.tsv"
    assert summary_path.exists()
    with open(summary_path) as f:
        lines = f.readlines()
        assert len(lines) == 4

    assert (outdir / "run.json").exists()
    # Check the number of transcripts in the HTML report
    html_path = outdir / "qc_report.html"
    assert html_path.exists()
    html = html_path.read_text()

    assert 'class="stats"' in html
    assert 'class="stat-box"' in html
    assert "Total Transcripts" in html
    assert ">3<" in html


if __name__ == "__main__":
    pytest.main()
