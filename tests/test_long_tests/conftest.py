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
    def _get_results_func(n_jobs):
        if True:
            print("PROCESSING!")
            primers = make_linear(load_fasta_glob(paths["primers"]))
            templates = load_genbank_glob(paths["templates"])

            query_path = join(here, "data/test_data/genbank/designs/*.gb")
            queries = make_circular(load_genbank_glob(query_path))[:LIM_NUM_DESIGNS]

            design = Design(span_cost=cached_span_cost)
            design.add_materials(primers=primers, templates=templates, queries=queries)
            if n_jobs > 1:
                design._run_with_pool(n_jobs, 1)
            else:
                design.run()
            return design, design.results

    return _get_results_func


@pytest.fixture(scope="session")
def single_processed_results(
    _processed_results
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _processed_results(1)


@pytest.fixture(scope="session")
def multi_processed_results(
    _processed_results
) -> Tuple[Design, Dict[str, DesignResult]]:
    return _processed_results(10)
