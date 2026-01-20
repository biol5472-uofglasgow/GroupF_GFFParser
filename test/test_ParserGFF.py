import pytest

from gene_summariser.praser import ParserGFF

@pytest.fixture
def parser():
    gff_path = "fixture_data/models.gff3"
    return ParserGFF(gff_path)

def test_parse_transcript(parser):
    transcripts = parser.parse_transcripts
    assert len(transcripts) > 0 