import pytest
from Bio import SeqIO

from gene_summariser.fasta import get_full_sequence
from gene_summariser.models import CDS, Transcript
from gene_summariser.parser import ParserGFF


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

    assert full_cds == ("TTGACGTACGT" "CGTACGTACGTACGTACGTAC" "TACGTACGTACGTTGC")


def test_full_transcript(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]

    genome = {
        record.id: record
        for record in SeqIO.parse("test/fixtures/testfasta.fasta", "fasta")
    }

    full_sequence = get_full_sequence(genome, transcript)
    assert len(full_sequence) % 3 == 0


def test_reverse_complement(parser):
    negative_transcript = Transcript(
        transcript_id="neg_tx",
        gene_id="neg_gene",
        seqid="chr1",
        start=1,
        end=10,
        strand="-",
        exons=[],
        cds_features=[
            CDS(seqid="chr1", start=1, end=10, strand="-"),
        ],
        attributes={},
    )

    genome = {
        record.id: record
        for record in SeqIO.parse("test/fixtures/testfasta.fasta", "fasta")
    }
    full_sequence = get_full_sequence(genome, negative_transcript)
    assert full_sequence == "ATACGTACGT"
