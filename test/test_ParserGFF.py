from gene_summariser.praser import ParserGFF
import pytest

@pytest.fixture
def gff_file():
    return "test/models.gff3"

@pytest.fixture
def parser(gff_file):
    return ParserGFF(gff_file)

def test_number_transcripts(parser):
    transcripts = parser.parse_transcripts()
    assert len(transcripts) == 2

def test_correct_feilds(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[1]
    
    assert transcript.transcript_id == "tx2"
    assert transcript.gene_id == "gene1"
    assert transcript.seqid == "chr1"
    assert transcript.start == 5
    assert transcript.end == 80
    assert transcript.strand == "+"
    assert len(transcript.exons) == 2
    assert len(transcript.cds_features) == 0
