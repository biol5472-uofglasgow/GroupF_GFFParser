from typing import Literal

import pytest

from gene_summariser.models import CDS, Transcript
from gene_summariser.parser import ParserGFF
from gene_summariser.qc import QCChecker


class TestFeatureGFF:
    def __init__(self, attributes, feature_id):
        self.attributes = attributes
        self.id = feature_id


@pytest.fixture
def gff_file():
    return "test/models.gff3"


@pytest.fixture
def parser(gff_file: Literal["test/models.gff3"]):
    return ParserGFF(gff_file)


def test_divisable_by_three(parser: ParserGFF):
    transcripts = parser.parse_transcripts()
    QC_checker = QCChecker("test/testfasta.fasta")

    assert "CDS_NOT_DEVISABLE_BY_3" not in QC_checker.check_transcript(transcripts[0])
    assert "CDS_NOT_DEVISABLE_BY_3" not in QC_checker.check_transcript(transcripts[1])


def test_divisable_not_divisable_by_three():
    QC_checker = QCChecker("test/testfasta.fasta")

    cds_features_fail = [
        CDS(seqid="chr1", start=1, end=10, strand="+"),
        CDS(seqid="chr1", start=11, end=14, strand="+"),
    ]
    transcript = Transcript(
        transcript_id="fail_test",
        gene_id="ID1",
        seqid="chr1",
        start=1,
        end=100,
        strand="+",
        exons=[],
        cds_features=cds_features_fail,
        attributes={},
    )

    assert "cds_not_divisible_by_3" in QC_checker.check_transcript(transcript)


def test_no_start_codon_flag(parser: ParserGFF):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "no_start_codon" not in flags


def test_stop_codon_flags(parser: ParserGFF):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "missing_stop_codon" in flags
