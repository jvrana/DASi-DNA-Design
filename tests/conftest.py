import functools
import hashlib
import os
import random
import warnings
from glob import glob
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Dict

import numpy as np
import pandas as pd
import pylab as plt
import pytest
from Bio import BiopythonParserWarning
from pyblast import BioBlastFactory
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import cost
from dasi.log import logger

##############################
# Global setup
##############################
warnings.simplefilter("ignore", BiopythonParserWarning)
logger.set_level("DEBUG")


def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line("markers", "slowtest: marks a long running test")
    config.addinivalue_line("markers", "benchmark: marks a benchmarking test")


def _pytest_auto_mark_benchmark(item):
    """Automatically mark tests that use the `benchmark` fixture."""
    marks = [mark for mark in item.iter_markers(name="benchmark")]
    if not marks and "benchmark" in item.funcargnames:
        item.add_marker(pytest.mark.benchmark)


def pytest_runtest_setup(item):
    _pytest_auto_mark_benchmark(item)


##############################
# Setup pandas display options
##############################

desired_width = 500
pd.set_option("display.width", desired_width)
# np.set_printoption(linewidth=desired_width)
pd.set_option("display.max_columns", 20)

##############################
# Test paths
##############################


@pytest.fixture(scope="session")
def here():
    return dirname(abspath(__file__))


@pytest.fixture(scope="session")
def fixture_dir(here):
    return join(here, "fixtures")


@pytest.fixture(scope="session")
def cost_filepath(fixture_dir):
    return join(fixture_dir, "span_cost.b")


@pytest.fixture(scope="session")
def cost_checksum_filepath(fixture_dir):
    return join(fixture_dir, "cost.checksum")


PRIMERS = "primers"
TEMPLATES = "templates"
QUERIES = "queries"
REGISTRY = "registry"
TEST_DATA = "test_data"


@pytest.fixture(scope="session")
def paths(here) -> Dict[str, str]:
    return {
        PRIMERS: join(here, "data/test_data/primers/primers.fasta"),
        TEMPLATES: join(here, "data/test_data/genbank/templates/*.gb"),
        QUERIES: join(here, "data/test_data/genbank/designs/*.gb"),
        REGISTRY: join(here, "data/test_data/genbank/benchling_registry/*.gb"),
        TEST_DATA: join(here, "data/test_data"),
    }


##############################
# Fixtures
##############################
@pytest.fixture(scope="session")
def blast_factory(paths) -> BioBlastFactory:
    factory = BioBlastFactory()

    primers = make_linear(load_fasta_glob(paths[PRIMERS]))
    templates = load_genbank_glob(paths[REGISTRY])
    queries = make_circular(load_genbank_glob(paths[QUERIES]))

    factory.add_records(primers, PRIMERS)
    factory.add_records(templates, TEMPLATES)
    factory.add_records(queries, QUERIES)

    return factory


##############################
# Cost Fixtures
##############################


def hashfiles(files, hash="sha1", hashfunc=None):
    if not hashfunc:

        def hashfunc(string):
            return getattr(hashlib, hash)(string.encode("utf-8")).hexdigest()

    contents = ""
    sorted_path = sorted(files)
    for file in sorted_path:
        with open(file) as f:
            contents += hashfunc(f.read())
    return hashfunc(contents)


def cost_module_checksum():
    cost_dir = os.path.dirname(cost.__file__)
    cost_files = sorted(
        glob(os.path.join(cost_dir, "*.py")) + glob(os.path.join(cost_dir, "*.json"))
    )
    return hashfiles(cost_files)


def cached(path, save_func, load_func, checksum_path, logger=None):
    different_checksum = True
    checksum = cost_module_checksum()
    if os.path.isfile(checksum_path):
        with open(checksum_path) as f:
            stored_checksum = f.read().strip()
            if stored_checksum == checksum:
                different_checksum = False
            if logger:
                logger.debug("Stored checksum: {}".format(stored_checksum))
    if logger:
        logger.debug("Checksum: {}".format(checksum))
    if different_checksum or not os.path.isfile(path):
        if logger:
            logger.debug("Using default params")
        model = save_func(path)
        with open(checksum_path, "w") as f:
            f.write(checksum)
        stat = os.stat(path)
        if logger:
            logger.debug("Wrote {} bytes".format(stat.st_size))
    else:
        if logger:
            logger.debug("Loading {}".format(path))
        stat = os.stat(path)
        if logger:
            logger.debug("Loaded {} bytes".format(stat.st_size))
        model = load_func(path)
    return model


@pytest.fixture(scope="session")
def cached_span_cost(cost_filepath, cost_checksum_filepath):
    """This will check the checksum of the cost module against the last
    checksum. If checksums are the same, the span cost will be loaded. Else,
    span_cost will be created from default parameters and saved with the cost
    module's checksum.

    :param cost_filepath: path of the span_cost
    :param cost_checksum_filepath: path of the checksum
    :return: SpanCost
    """

    def load_span_cost(path):
        span_cost = cost.SpanCost.load(path)
        return span_cost

    def save_span_cost(path):
        span_cost = cost.SpanCost.open()
        span_cost.dump(path)
        return span_cost

    return cached(
        cost_filepath,
        load_func=load_span_cost,
        save_func=save_span_cost,
        checksum_path=cost_checksum_filepath,
        logger=logger,
    )


span_cost = cached_span_cost

##############################
# Auto Fixtures
##############################


@pytest.fixture(autouse=True, scope="function")
def random_seed(request):
    seed = request.param_index
    print("RANDOM SEED({})".format(seed))
    random.seed(seed)
    np.random.seed(seed)


@pytest.fixture(autouse=True)
def patch_matplotlib_show(tmpdir, monkeypatch):
    def patched_show(*args, **kwargs):
        path = tmpdir.join("image.png")
        plt.savefig(path, format="png")

    monkeypatch.setattr(plt, "show", patched_show)
