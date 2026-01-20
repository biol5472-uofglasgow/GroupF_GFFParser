import pytest

from gene_summariser.praser import ParserGFF

@pytest.fixture
def parser():
    gff_path = "fixture_data/models.gff3"
    return ParserGFF(gff_path)