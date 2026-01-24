import pytest

from gene_summariser.parser import ParserGFF


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


def test_exon_parsing_gene_id(parser):
    feature = TestFeatureGFF(attributes={"gene_id": "gene1"}, feature_id="tx1")
    gene_id = ParserGFF._get_id(feature)
    assert gene_id == "gene1"


def test_exon_parsing_parent(parser):
    feature = TestFeatureGFF(attributes={"Parent": "gene2"}, feature_id="tx2")
    gene_id = ParserGFF._get_id(feature)
    assert gene_id == "gene2"


def test_exon_parsing_fallback(parser):
    feature = TestFeatureGFF(attributes={}, feature_id="gene3")
    gene_id = ParserGFF._get_id(feature)
    assert gene_id == "gene3"


def test_cds_parsing(parser):
    transcripts = parser.parse_transcripts()
    transcript = transcripts[0]

    assert len(transcript.cds_features) == 3
    assert transcript.cds_features[0].start == 10
    assert transcript.cds_features[0].end == 20
    assert transcript.cds_features[0].strand == "+"
    assert transcript.cds_features[0].phase == 0

    assert transcript.cds_features[1].start == 30
    assert transcript.cds_features[1].end == 50
    assert transcript.cds_features[1].strand == "+"
    assert transcript.cds_features[1].phase == 2

    assert transcript.cds_features[2].start == 60
    assert transcript.cds_features[2].end == 75
    assert transcript.cds_features[2].strand == "+"
    assert transcript.cds_features[2].phase == 1
