from dasi.design import Design, DesignResult
from dasi.alignments import AlignmentContainer
from dasi.cost import SpanCost
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
import pytest
from os.path import join
from more_itertools import pairwise
import networkx as nx
import pickle


@pytest.mark.parametrize(
    "query", ["pmodkan-ho-pact1-z4-er-vpr.gb", "plko-pef1a-frt-tdtomato-wpre.gb"]
)
def test_num_groups_vs_endpoints(here, paths, query, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design._blast()
    containers = design.container_list()
    assert len(containers) == 1
    container = containers[0]
    container.expand()
    groups = container.groups()
    print(len(groups) ** 2)

    a_arr = set()
    b_arr = set()

    for g in groups:
        a_arr.add(g.query_region.a)
        b_arr.add(g.query_region.b)

    print(len(a_arr) * len(b_arr))


def print_edge_cost(path, graph):
    total = 0
    path = path[:] + path[:1]
    for n1, n2 in pairwise(path):
        try:
            edata = graph[n1][n2]
            total += edata["weight"]
            print((n1, n2, edata["weight"]))
        except:
            print((n1, n2, "MISSING EDGE"))

    print("TOTAL: {}".format(total))


@pytest.mark.parametrize(
    "query",
    [
        "pmodkan-ho-pact1-z4-er-vpr.gb",
        "plko-pef1a-frt-tdtomato-wpre.gb",
        "pins-01-hu6-sv40-nt1-optgrna.gb",
        "pins-01-hu6-r1-optgrna.gb",
        "pins-01-hu6-r5-optgrna.gb",
    ],
)
def test_real_design(here, paths, query, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) == 1

    results = design.optimize()

    for query_key, result in results.items():
        assembly = result.assemblies[0]
        assembly.print()

        # df = assembly.to_df()
        # filename = "{}__{}__test__.csv".format(assembly.query_key, assembly.query.name)
        # df.to_csv(filename)

        assert len(result.query) == sum(assembly.to_df()["span"])


@pytest.mark.parametrize("query", ["pmodkan-ho-pact1-z4-er-vpr.gb"])
def test_profile_compile(here, paths, query, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)
    design.compile()


def test_real_design2(here, paths, span_cost):
    query = "goal1.gb"
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["registry"])

    fragments = [
        f for f in templates if f.annotations.get("topology", None) == "linear"
    ]
    plasmids = [
        f for f in templates if f.annotations.get("topology", None) == "circular"
    ]

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_fragments(fragments)
    design.add_templates(fragments + plasmids)
    design.add_primers(primers)
    design.add_queries(queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) == 1

    results = design.optimize()

    for query_key, result in results.items():
        assembly = result.assemblies[0]

        assert len(result.query) == sum(assembly.to_df()["span"])

        assembly.print()


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


def test_multidesign(here, paths, span_cost):
    """Expect more than one graph to be output if multiple queries are provided"""
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) > 1

    design.optimize()


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
            for c in df.columns:
                assert np.any(df[c].to_numpy().flatten())
            # assert df['subject'].to_numpy().flatten()
            # assert df['subject_ends'].to_numpy().flatten()
            print(assembly.to_df())


def test_pickle_design_result():
    container = AlignmentContainer({'none', None}, [])
    pickle.dumps(container)
    result = DesignResult(AlignmentContainer({'none': None}, []), nx.DiGraph, "none")
    s = pickle.dumps(result)
    unpickled_result = pickle.loads(s)