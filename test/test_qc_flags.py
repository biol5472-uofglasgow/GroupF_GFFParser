from typing import Literal

from gene_summariser.models import CDS, Transcript
import pytest

from gene_summariser.models import CDS, Transcript
from gene_summariser.parser import ParserGFF
from gene_summariser.qc import QCChecker


class TestFeatureGFF:
    def __init__(self, attributes, feature_id):
        self.attributes = attributes
        self.id = feature_id


# Set up fixtures for GFF file and parser
@pytest.fixture
def gff_file():
    return "test/models.gff3"


@pytest.fixture
def parser(gff_file: Literal["test/models.gff3"]):
    return ParserGFF(gff_file)


def test_divisable_by_three(parser: ParserGFF):
    """
    Unit test for the following condition:
    - A transcript with CDS features that are all divisible by 3 does not have the "cds_not_divisible_by_3" flag.
    - A transcript with no CDS features does not have the "cds_not_divisible_by_3" flag.
    """
    transcripts = parser.parse_transcripts()
    QC_checker = QCChecker("test/testfasta.fasta")

    assert "cds_not_divisible_by_3" not in QC_checker.check_transcript(transcripts[0])
    assert "cds_not_divisible_by_3" not in QC_checker.check_transcript(transcripts[1])


def test_divisable_not_divisable_by_three():
    """
    Unit test for the following condition:
    - A transcript with CDS features that are not all divisible by 3 has the "cds_not_divisible_by_3" flag.
    """
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
    """
    Unit test for the following condition:
    - A transcript with start codons present does not have the "no_start_codon" flag.
    """
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "no_start_codon" not in flags


def test_stop_codon_flags(parser: ParserGFF):
    """
    Unit test for the following condition:
    - A transcript missing a stop codon has the "missing_stop_codon" flag.
    """
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "missing_stop_codon" in flags


def test_cds_phase_consistent_pass():
    """
    Transcript with correct CDS phase progression
    should NOT be flagged.
    """
    qc_checker = QCChecker(fasta_file=None)
    
    correct_cds= [
        CDS(seqid="chr1", start=1, end=100, strand="+", phase=0),
        CDS(seqid="chr1", start=201, end=300, strand="+", phase=1),  # CORRECT
    ]

    correct_transcript = Transcript(
        transcript_id="tx_pass",
        gene_id="gene1",
        seqid="chr1",
        start=1,
        end=300,
        strand="+",
        exons=[],
        cds_features=correct_cds,
        attributes={},
    )

    
    flags_correct = qc_checker.check_transcript(correct_transcript)
    assert "cds_phase_inconsistent" not in flags_correct

    wrong_cds = [
        CDS(seqid="chr1", start=1, end=100, strand="+", phase=0),
        CDS(seqid="chr1", start=201, end=300, strand="+", phase=0),  # wrong
    ]
    wrong_transcript = Transcript(
        transcript_id="tx_bad",
        gene_id="gene1",
        seqid="chr1",
        start=1,
        end=300,
        strand="+",
        exons=[],
        cds_features=wrong_cds,
        attributes={},
    )

    flags_bad = qc_checker.check_transcript(wrong_transcript)
    assert "cds_phase_inconsistent" in flags_bad
