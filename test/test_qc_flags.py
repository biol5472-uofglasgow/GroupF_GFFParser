import pytest

from gene_summariser.models import CDS, Transcript
from gene_summariser.parser import ParserGFF
from gene_summariser.qc import append_flags_to_summary


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


def test_divisable_by_three(parser):
    transcripts = parser.parse_transcripts()
    summaries = append_flags_to_summary(transcripts)

    assert "CDS_NOT_DEVISABLE_BY_3" not in summaries[0].flags
    assert "CDS_NOT_DEVISABLE_BY_3" not in summaries[1].flags


def test_divisable_not_divisable_by_three():
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

    transcripts = [transcript]
    summaries = append_flags_to_summary(transcripts)
    assert "CDS_NOT_DEVISABLE_BY_3" in summaries[0].flags
