import pytest

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
def parser(gff_file):
    return ParserGFF(gff_file)


# def test_divisable_by_three(parser):
#     transcripts = parser.parse_transcripts()
#     QC_checker = QCChecker('test/')

#     assert "CDS_NOT_DEVISABLE_BY_3" not in summaries[0].flags
#     assert "CDS_NOT_DEVISABLE_BY_3" not in summaries[1].flags


# def test_divisable_not_divisable_by_three():
#     cds_features_fail = [
#         CDS(seqid="chr1", start=1, end=10, strand="+"),
#         CDS(seqid="chr1", start=11, end=14, strand="+"),
#     ]
#     transcript = Transcript(
#         transcript_id="fail_test",
#         gene_id="ID1",
#         seqid="chr1",
#         start=1,
#         end=100,
#         strand="+",
#         exons=[],
#         cds_features=cds_features_fail,
#         attributes={},
#     )

#     transcripts = [transcript]
#     summaries = append_flags_to_summary(transcripts)
#     assert "CDS_NOT_DEVISABLE_BY_3" in summaries[0].flags


def test_no_start_codon_flag(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "no_start_codon" not in flags


def test_stop_codon_flags(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]  # Transcript with no CDS features
    QC_checker = QCChecker("test/testfasta.fasta")
    flags = QC_checker.check_transcript(transcript)
    assert "missing_stop_codon" in flags
