import pytest
from dasi import Design
from dasi.cost import SpanCost
from dasi.utils import multipoint_shortest_path
from pyblast.utils import make_linear, make_circular, load_genbank_glob, load_fasta_glob
from os.path import join, isfile, dirname, abspath
import networkx as nx
import numpy as np
import dill
from dasi.log import logger
from more_itertools import pairwise


@pytest.fixture(scope='module')
def span_cost():
    here = dirname(abspath(__file__))
    filepath = join(here, 'span_cost.pkl')
    if isfile(filepath):
        with open(filepath, 'rb') as f:
            span_cost = dill.load(f)
    else:
        span_cost = SpanCost()
        with open(filepath, 'wb') as f:
            dill.dump(span_cost, f)
    return span_cost


@pytest.mark.parametrize('query', [
    "pmodkan-ho-pact1-z4-er-vpr.gb",
    'plko-pef1a-frt-tdtomato-wpre.gb',
    'pins-01-hu6-sv40-nt1-optgrna.gb'
])
def test_3point_optimization(here, paths, query, span_cost):
    with logger.timeit("INFO", "Loading files"):
        primers = make_linear(load_fasta_glob(paths["primers"]))
        templates = load_genbank_glob(paths["templates"])
        query_path = join(here, 'data/test_data/genbank/designs', query)
        queries = make_circular(load_genbank_glob(query_path))

    with logger.timeit("INFO", "Compiling"):
        design = Design(span_cost)
        design.add_materials(primers=primers, templates=templates, queries=queries)
        design.compile()

    graph = list(design.graphs.values())[0]

    nodelist = list(graph.nodes())
    node_to_i = {v: i for i, v in enumerate(nodelist)}

    with logger.timeit("INFO", "FloydWarshall"):
        weight_matrix = np.array(nx.floyd_warshall_numpy(graph, nodelist=nodelist, weight='weight'))

    # collect all 3-point cycles
    with logger.timeit("INFO", "3PointOptimization") as log:
        cycles = []
        for i, n1 in enumerate(nodelist):
            if n1[2] != 'A':
                continue
            for n2 in graph.successors(n1):
                j = node_to_i[n2]
                if i == j:
                    continue
                a = weight_matrix[i, j]
                if a == np.inf:
                    continue
                for k in range(len(weight_matrix[0])):
                    if k == i or k == j:
                        continue
                    n3 = nodelist[k]
                    if n3[2] != 'A':
                        continue
                    b = weight_matrix[j, k]
                    if b == np.inf:
                        continue

                    c = weight_matrix[k, i]
                    if c == np.inf:
                        continue

                    if a + b + c != np.inf:
                        x = ((n1, n2, n3), (a, b, c), a + b + c)
                        cycles.append(x)
        log.info("Found {} cycles".format(len(cycles)))
        cycles = sorted(cycles, key=lambda x: x[-1])

    with logger.timeit("INFO", "CollectingPaths") as log:
        n_paths = 20
        unique_cyclic_paths = []
        for c in cycles:
            if len(unique_cyclic_paths) >= n_paths:
                break
            path = multipoint_shortest_path(graph, c[0], weight_key='weight', cyclic=True)
            if path not in unique_cyclic_paths:
                unique_cyclic_paths.append(path)
        print(unique_cyclic_paths)

    best_path = unique_cyclic_paths[0]

    container = design.container_list()[0]

    for n1, n2 in pairwise(best_path + best_path[:1]):
        a, b = n1[0], n2[0]
        edata = graph[n1][n2]
        print(edata['weight'])
        groups = container.find_group_by_pos(a, b)
        pass