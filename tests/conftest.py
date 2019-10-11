import warnings
from os.path import abspath
from os.path import dirname
from os.path import join

import pytest
from Bio import BiopythonParserWarning
from pyblast import BioBlastFactory
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.log import logger

warnings.simplefilter("ignore", BiopythonParserWarning)

logger.set_level("DEBUG")


import pandas as pd

desired_width = 500

pd.set_option("display.width", desired_width)

# np.set_printoption(linewidth=desired_width)

pd.set_option("display.max_columns", 20)


@pytest.fixture(scope="session")
def here():
    return dirname(abspath(__file__))


PRIMERS = "primers"
TEMPLATES = "templates"
QUERIES = "queries"
REGISTRY = "registry"


@pytest.fixture(scope="session")
def paths(here):
    return {
        PRIMERS: join(here, "data/test_data/primers/primers.fasta"),
        TEMPLATES: join(here, "data/test_data/genbank/templates/*.gb"),
        QUERIES: join(
            here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb"
        ),
        REGISTRY: join(here, "data/test_data/genbank/benchling_registry/*.gb"),
    }


@pytest.fixture(scope="session")
def blast_factory(paths):
    factory = BioBlastFactory()

    primers = make_linear(load_fasta_glob(paths[PRIMERS]))
    templates = load_genbank_glob(paths[TEMPLATES])
    queries = make_circular(load_genbank_glob(paths[QUERIES]))

    factory.add_records(primers, PRIMERS)
    factory.add_records(templates, TEMPLATES)
    factory.add_records(queries, QUERIES)

    return factory


def pytest_assertrepr_compare(config, op, left, right):
    if op in ("==", "!="):
        return ["{} {} {}".format(left, op, right)]
