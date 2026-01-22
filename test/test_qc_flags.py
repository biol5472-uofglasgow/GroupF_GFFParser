from gene_summariser.praser import ParserGFF
from gene_summariser.qc import append_flags_to_summary, flag_devisable_by_three
import pytest

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

def test_devisable_by_three(parser):
    transcripts = parser.parse_transcripts()
    summary1 = append_flags_to_summary([transcripts[0]])
    assert "CDS_NOT_DEVISABLE_BY_3" not in summary1.flags

    summary2 = append_flags_to_summary([transcripts[1]])
    assert "CDS_NOT_DEVISABLE_BY_3" in summary2.flags