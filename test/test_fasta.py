import pytest

from gene_summariser.fasta import iter_cds_sequences
from gene_summariser.parser import ParserGFF


@pytest.fixture
def gff_file():
    return "test/models.gff3"


@pytest.fixture
def parser(gff_file):
    return ParserGFF(gff_file)


def test_extract_gene_sequences(parser):
    transcripts = parser.parse_transcripts()
    fasta_file = "test/testfasta.fasta"
    transcript = transcripts[0]
    cds_parts = list(iter_cds_sequences(fasta_file, transcript))

    assert cds_parts[0] == "CGTACGTACGT"
    assert cds_parts[1] == "TACGTACGTACGTACGTAC"
    assert cds_parts[2] == "ACGTACGTACGTACG"


def test_phase(parser):
    transcripts = parser.parse_transcripts()
    fasta_file = "test/testfasta.fasta"
    transcript = transcripts[0]
    cds_parts = list(iter_cds_sequences(fasta_file, transcript))

    assert len(cds_parts[1]) == (50 - 30 + 1) - 2


def test_full_transcript(parser):
    transcripts = parser.parse_transcripts()
    fasta_file = "test/testfasta.fasta"
    transcript = transcripts[0]
    cds_parts = list(iter_cds_sequences(fasta_file, transcript))

    full_sequence = "".join(cds_parts)
    assert len(full_sequence) % 3 == 0
