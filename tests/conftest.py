import pytest
from os.path import dirname, abspath, join
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
from pyblast import BioBlastFactory
from dasi.log import logger
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

logger.set_level("INFO")


@pytest.fixture(scope="module")
def here():
    return dirname(abspath(__file__))


PRIMERS = "primers"
TEMPLATES = "templates"
QUERIES = "queries"
REGISTRY = "registry"


@pytest.fixture(scope="module")
def paths(here):
    return {
        PRIMERS: join(here, "data/test_data/primers/primers.fasta"),
        TEMPLATES: join(here, "data/test_data/genbank/templates/*.gb"),
        QUERIES: join(
            here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb"
        ),
        REGISTRY: join(here, 'data/test_data/genbank/benchling_registry/*.gb')
    }


@pytest.fixture(scope="module")
def blast_factory(paths):
    factory = BioBlastFactory()

    primers = make_linear(load_fasta_glob(paths[PRIMERS]))
    templates = load_genbank_glob(paths[TEMPLATES])
    queries = make_circular(load_genbank_glob(paths[QUERIES]))

    factory.add_records(primers, PRIMERS)
    factory.add_records(templates, TEMPLATES)
    factory.add_records(queries, QUERIES)

    return factory
