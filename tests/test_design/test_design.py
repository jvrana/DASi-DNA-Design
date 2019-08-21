from dasi.design import Design
from dasi.cost import SpanCost
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
import pytest
from os.path import join
from more_itertools import pairwise


span_cost = SpanCost()


@pytest.mark.parametrize('query', [
    "pmodkan-ho-pact1-z4-er-vpr.gb",
    'plko-pef1a-frt-tdtomato-wpre.gb'
])
def test_num_groups_vs_endpoints(here, paths, query):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, 'data/test_data/genbank/designs', query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design._blast()
    containers = design.container_list()
    assert len(containers) == 1
    container = containers[0]
    container.expand()
    groups = container.groups()
    print(len(groups)**2)

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
            total += edata['weight']
            print((n1, n2, edata['weight']))
        except:
            print((n1, n2, "MISSING EDGE"))

    print("TOTAL: {}".format(total))

@pytest.mark.parametrize('query', [
    # "pmodkan-ho-pact1-z4-er-vpr.gb",
    'plko-pef1a-frt-tdtomato-wpre.gb',
    # 'pins-01-hu6-sv40-nt1-optgrna.gb'
])
def test_real_design(here, paths, query):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, 'data/test_data/genbank/designs', query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) == 1
    path_dict = design.optimize()
    best_path = list(path_dict.values())[0][0]
    df = design.design()
    print(df)

    print_edge_cost(best_path, list(design.graphs.values())[0])


def test_multidesign(here, paths):
    """Expect more than one graph to be output if multiple queries are provided"""
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, 'data/test_data/genbank/designs/*.gb')
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) > 1

    design.optimize()
