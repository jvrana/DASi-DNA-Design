import os
from copy import deepcopy
from os.path import join

import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.cost import SpanCost
from dasi.design import Design
from dasi.log import logger


here = os.path.abspath(os.path.dirname(__file__))
do_save = True


@pytest.fixture(scope="session")
def span_cost():
    """Saves the span cost as bytes; reloads when called."""
    path = os.path.join(here, "span_cost.b")
    if do_save and os.path.isfile(path):
        with logger.timeit("INFO", "loading bytes"):
            print("Loading file: {}".format(path))
            span_cost = SpanCost.load(path)
    else:
        span_cost = SpanCost.default()
        if do_save:
            with logger.timeit("INFO", "saving bytes"):
                print("Saving file: {}".format(path))
                span_cost.dump(path)
    return span_cost


@pytest.fixture(scope="session")
def _processed_results(here, paths, span_cost):
    def _get_results_func(n_jobs, only_compile=False, precompiled=None):
        if precompiled:
            design = precompiled
        else:
            print("PROCESSING!")
            primers = make_linear(load_fasta_glob(paths["primers"]))
            templates = load_genbank_glob(paths["templates"])

            query_path = join(here, "data/test_data/genbank/designs/*.gb")
            queries = make_circular(load_genbank_glob(query_path))[:3]

            design = Design(span_cost=span_cost)
            design.add_materials(primers=primers, templates=templates, queries=queries)
            design.n_jobs = n_jobs
            design.compile()
        if only_compile:
            return design, {}
        else:
            print("OPTIMIZING!")
            for container in design.containers.values():
                print(len(container.groups()))
            results = design.optimize()
            return design, results

    return _get_results_func


@pytest.fixture(scope="session")
def _single_compiled_design(_processed_results):
    return _processed_results(1, only_compile=True)[0]


@pytest.fixture(scope="session")
def _multi_compiled_design(_processed_results):
    return _processed_results(10, only_compile=True)[0]


@pytest.fixture(scope="session")
def _single_processed_results(_processed_results, _single_compiled_design):
    return _processed_results(1, precompiled=_single_compiled_design)


@pytest.fixture(scope="session")
def _multi_processed_results(_processed_results, _multi_compiled_design):
    return _processed_results(10, precompiled=_multi_compiled_design)


@pytest.fixture(scope="session")
def single_compiled_results(_single_compiled_design):
    return deepcopy(_single_compiled_design)


@pytest.fixture(scope="session")
def multi_compiled_results(_multi_compiled_design):
    return deepcopy(_multi_compiled_design)


@pytest.fixture(scope="session")
def single_processed_results(_single_processed_results):
    return deepcopy(_single_processed_results)


@pytest.fixture(scope="session")
def multi_processed_results(_multi_processed_results):
    return deepcopy(_multi_processed_results)
