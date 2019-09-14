import pytest
from pyblast.utils import make_linear, make_circular, load_fasta_glob, load_genbank_glob
from os.path import join
from dasi import Design
from dasi.design import DesignResult
from dasi.alignments import AlignmentContainer
import pickle
import networkx as nx


def test_pickle_design_result():
    container = AlignmentContainer({'none', None}, [])
    pickle.dumps(container)
    result = DesignResult(AlignmentContainer({'none': None}, []), nx.DiGraph, "none")
    s = pickle.dumps(result)
    unpickled_result = pickle.loads(s)

def f(arg):
    scost, primers, templates, queries, results = arg
    design = Design(span_cost=scost)
    design.add_materials(primers=primers, templates=templates, queries=queries)
    design.compile()
    return design.optimize()

@pytest.mark.parametrize("ncores", [10])
def test_multiprocessing(here, paths, span_cost, ncores):
    """Test that demonstrates how multiprocessing can speed up designing multiple constructs."""
    from multiprocessing import Pool

    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))

    args = [(span_cost, primers, templates, [query], []) for query in queries]
    print("Number of queries: {}".format(len(queries)))
    with Pool(processes=ncores) as pool:  # start 4 worker processes
        results = pool.map(f, args)
    print(results)
    assert results

# # long tests
# def test_multidesign(here, paths, span_cost):
#     """Expect more than one graph to be output if multiple queries are provided"""
#     primers = make_linear(load_fasta_glob(paths["primers"]))
#     templates = load_genbank_glob(paths["templates"])
#
#     query_path = join(here, "data/test_data/genbank/designs/*.gb")
#     queries = make_circular(load_genbank_glob(query_path))
#
#     design = Design(span_cost=span_cost)
#
#     design.add_materials(primers=primers, templates=templates, queries=queries)
#
#     design.compile()
#
#     assert len(design.graphs) == len(queries)
#     assert len(design.graphs) > 1
#
#     design.optimize()


def test_multiprocessing_multidesign(here, paths, span_cost):
    """Expect more than one graph to be output if multiple queries are provided"""
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)
    design.n_jobs = 10
    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) > 1

    results = design.optimize()

    # assert
    import numpy as np
    for result in results.values():
        for assembly in result.assemblies:
            df = assembly.to_df()
            print(assembly.to_df())
            for c in df.columns:
                assert np.any(df[c].to_numpy().flatten())
            # assert df['subject'].to_numpy().flatten()
            # assert df['subject_ends'].to_numpy().flatten()
