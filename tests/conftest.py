import pytest
from os.path import dirname, abspath, join
from pyblast import BlastBase, BioBlast
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear
from shoestring import BioBlastFactory


@pytest.fixture(scope="module")
def here():
    return dirname(abspath(__file__))

PRIMERS = "primers"
TEMPLATES = "templates"
QUERIES = "queries"

@pytest.fixture(scope="module")
def blast_factory(here):
    factory = BioBlastFactory()

    primers = make_linear(
        load_fasta_glob(join(here, "data/test_data/primers/primers.fasta"))
    )
    templates = load_genbank_glob(join(here, "data/test_data/genbank/templates/*.gb"))
    queries = load_genbank_glob(
        join(here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb")
    )

    factory.add_records(primers, PRIMERS)
    factory.add_records(templates, TEMPLATES)
    factory.add_records(queries, QUERIES)

    return (
        factory
    )
