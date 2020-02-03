"""test_pickling.py.

Tests whether DesignResults can be pickled for multiprocessing
"""
import pickle

import networkx as nx
import pytest

from dasi import Design
from dasi import LibraryDesign
from dasi.design import DesignResult
from dasi.log import logger
from dasi.models import AlignmentContainer


def test_pickle_design_result():
    container = AlignmentContainer({"none", None}, [])
    pickle.dumps(container)
    result = DesignResult(AlignmentContainer({"none": None}, []), nx.DiGraph, "none")
    s = pickle.dumps(result)
    unpickled_result = pickle.loads(s)
    assert unpickled_result


@pytest.mark.slow
def test_large_pkl(span_cost):
    """Expect more than one graph to be output if multiple queries are
    provided."""
    design = Design.fake(
        n_designs=3, n_cyclic_seqs=100, n_linear_seqs=100, n_primers=100
    )
    design.compile()

    with logger.timeit("DEBUG", "pickling graphs"):
        pickle.loads(pickle.dumps(design.graphs))

    with logger.timeit("DEBUG", "pickling containers"):
        pickle.loads(pickle.dumps(design.container_factory))

    with logger.timeit("DEBUG", "pickling span_cost"):
        pickle.loads(pickle.dumps(span_cost))
