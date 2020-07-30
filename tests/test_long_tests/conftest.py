import os
from copy import deepcopy
from os.path import join
from typing import Callable
from typing import Dict
from typing import Tuple

import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.design import Design
from dasi.design import DesignResult

here = os.path.abspath(os.path.dirname(__file__))
do_save = True


LIM_NUM_DESIGNS = None


@pytest.fixture(scope="session")
def _processed_results(here, paths, cached_span_cost) -> Callable:
    def _get_results_func(n_jobs, only_compile=False, precompiled=None):
        if precompiled:
            design = precompiled
        else:
            print("PROCESSING!")
            primers = make_linear(load_fasta_glob(paths["primers"]))
            templates = load_genbank_glob(paths["templates"])

            query_path = join(here, "data/test_data/genbank/designs/*.gb")
            queries = make_circular(load_genbank_glob(query_path))[:LIM_NUM_DESIGNS]

            design = Design(span_cost=cached_span_cost)
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
def _single_compiled_design(_processed_results) -> Design:
    return _processed_results(1, only_compile=True)[0]


@pytest.fixture(scope="session")
def _multi_compiled_design(_processed_results) -> Design:
    return _processed_results(10, only_compile=True)[0]


@pytest.fixture(scope="session")
def _single_processed_results(
    _processed_results, _single_compiled_design
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _processed_results(1, precompiled=_single_compiled_design)


@pytest.fixture(scope="session")
def _multi_processed_results(
    _processed_results, _multi_compiled_design
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _processed_results(10, precompiled=_multi_compiled_design)


@pytest.fixture(scope="session")
def single_compiled_results(_single_compiled_design) -> Design:
    return _single_compiled_design


@pytest.fixture(scope="session")
def multi_compiled_results(_multi_compiled_design) -> Design:
    return _multi_compiled_design


@pytest.fixture(scope="session")
def single_processed_results(
    _single_processed_results,
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _single_processed_results


@pytest.fixture(scope="session")
def multi_processed_results(
    _multi_processed_results,
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _multi_processed_results
