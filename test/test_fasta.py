import pytest

from Bio import SeqIO
from gene_summariser.parser import ParserGFF\
from gene_summariser.fasta import get_full_sequence


@pytest.fixture
def gff_file():
    return "test/fixtures/models.gff3"


@pytest.fixture
def parser(gff_file):
    return ParserGFF(gff_file)


def test_extract_gene_sequences(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]
    genome = {
        record.id: record
        for record in SeqIO.parse("test/fixtures/testfasta.fasta", "fasta")
    }

    full_cds = get_full_sequence(genome, transcript)

    assert full_cds == (
        "TTGACGTACGT"
        "TACGTACGTACGTACGTAC"
        "ACGTACGTACGTTGC"
    )


def test_phase(parser):
    transcripts = parser.parse_transcripts()
    fasta_file = "test/fixtures/testfasta.fasta"
    transcript = transcripts[0]
    cds_parts = list(iter_cds_sequences(fasta_file, transcript))

    assert len(cds_parts[1]) == (50 - 30 + 1) - 2


def test_full_transcript(parser):
    transcripts = parser.parse_transcripts()
    fasta_file = "test/fixtures/testfasta.fasta"
    transcript = transcripts[0]
    cds_parts = list(iter_cds_sequences(fasta_file, transcript))

    full_sequence = "".join(cds_parts)
    assert len(full_sequence) % 3 == 0
