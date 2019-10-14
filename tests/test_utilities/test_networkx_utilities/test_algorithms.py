import random

import networkx as nx
import numpy as np
import pytest

from dasi.utils.networkx import find_all_min_paths
from dasi.utils.networkx import floyd_warshall_with_efficiency
from dasi.utils.networkx import sympy_dijkstras
from dasi.utils.networkx import sympy_floyd_warshall
from dasi.utils.networkx import sympy_multisource_dijkstras


COMPARISON_THRESHOLD = (
    0.01
)  # within 1% due to floating point errors during accumulating of multiple floats
EDGE_COMPARISON_THRESHOLD = 0.05


def add_data(G, u, v, weight, eff, weight_key="weight", eff_key="eff"):
    G.get_edge_data(u, v)[weight_key] = weight
    G.get_edge_data(u, v)[eff_key] = eff


def simple_graph():
    G = nx.path_graph(5)

    add_data(G, 0, 1, 100, 0.5)
    add_data(G, 1, 2, 100, 0.5)
    add_data(G, 2, 3, 100, 0.5)
    add_data(G, 3, 4, 100, 0.5)

    G.add_node(5)
    G.add_node(6)
    G.add_node(7)
    G.add_edge(0, 5)
    G.add_edge(5, 6)
    G.add_edge(6, 7)
    G.add_edge(7, 4)

    add_data(G, 0, 5, 200, 0.75)
    add_data(G, 5, 6, 200, 0.75)
    add_data(G, 6, 7, 200, 0.75)
    add_data(G, 7, 4, 200, 0.75)

    nodelist = list(G.nodes)
    return G, nodelist


def add_edata(G, weight="weight", eff="eff"):
    for n1, n2, edata in G.edges(data=True):
        edata[weight] = float(random.randint(0, 100))
        edata[eff] = np.random.uniform(0.8, 0.9)


def complete_graph(n_nodes, create_using=None):
    G = nx.complete_graph(n_nodes, create_using=create_using)
    add_edata(G)
    nodelist = list(G.nodes())
    return G, nodelist


def grid_graph(dims):
    G = nx.grid_graph(dims)
    add_edata(G)
    nodelist = list(G.nodes())
    return G, nodelist


def compare_floats(expected, x):
    diff = abs(expected - x)
    if diff / expected > EDGE_COMPARISON_THRESHOLD:
        return False
    return True


def compare(G, C, nodelist):
    d = find_all_min_paths(G)
    errors = []
    total = 0
    true_total = 0
    for u, v in np.ndenumerate(C):
        true_cost = d[nodelist[u[0]]][nodelist[u[1]]]["cost"]
        msg = "{} -> {} {} {}".format(u[0], u[1], true_cost, v)
        true_total += true_cost
        total += v
        if true_cost != v:
            if compare_floats(true_cost, v):
                status = "WARNING"
            else:
                status = "ERROR (T={})".format(EDGE_COMPARISON_THRESHOLD)
                errors.append(msg)
        else:
            status = ""

        print(status + " " + msg)
    if true_total == 0:
        if total:
            raise Exception("Threshold too large")
    else:
        diff = abs(total - true_total) / true_total
        if diff > COMPARISON_THRESHOLD:
            raise Exception("Threshold too large: {}".format(diff))
    assert not errors


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


class TestAllPairShortestPath:
    def floyd_warshall_compare(self, G, nodelist):
        C = floyd_warshall_with_efficiency(G, "weight", "eff", nodelist=nodelist)
        compare(G, C, nodelist)
        return C

    def test_simple(self):
        self.floyd_warshall_compare(*simple_graph())

    @pytest.mark.parametrize("n", [8])
    @pytest.mark.parametrize("repeat", range(3))
    def test_complete(self, n, repeat):
        C = self.floyd_warshall_compare(*complete_graph(n))
        assert check_symmetric(C)

    @pytest.mark.parametrize("n", list(range(2, 8)))
    @pytest.mark.parametrize("repeat", range(3))
    def test_complete_directed(self, n, repeat):
        """Checks to make sure the return result is not symmetric."""
        G, nodelist = complete_graph(n, create_using=nx.DiGraph)
        C = self.floyd_warshall_compare(G, nodelist)
        compare(G, C, nodelist)
        assert not check_symmetric(C)

    @pytest.mark.parametrize("using", [nx.Graph, nx.DiGraph])
    @pytest.mark.parametrize("n", [8])
    def test_equivalent_to_floyd_marshall(self, n, using):
        """Check in the trivial case where eff=1., that the result is
        equivalent to the result of the floyd_warshall_numpy algorithm."""
        G, nodelist = complete_graph(n, create_using=using)
        for n1, n2, edata in G.edges(data=True):
            edata["eff"] = 1.0
        M1 = nx.floyd_warshall_numpy(G, nodelist=nodelist, weight="weight")
        M2 = floyd_warshall_with_efficiency(G, "weight", "eff", nodelist=nodelist)
        assert np.allclose(M1, M2)

    def test_complete_large(self):
        self.floyd_warshall_compare(*complete_graph(9))

    @pytest.mark.parametrize(
        "dims", [[1, 1, 1], [2, 2, 2], [2, 3, 2], [5, 1, 2], [4, 2, 2]]
    )
    def test_grid_graph(self, dims):
        self.floyd_warshall_compare(*grid_graph(dims))

    def test_very_large_graph(self):
        G, nodelist = complete_graph(200)
        floyd_warshall_with_efficiency(G, "weight", "eff", nodelist=nodelist)

    def test_different_keys(self):
        G = nx.path_graph(10)
        add_edata(G, "cost", "efficiency")

        C = floyd_warshall_with_efficiency(G, weight_key="y", eff_key="x")
        assert 2.0 in C

        D = floyd_warshall_with_efficiency(G, weight_key="cost", eff_key="x")
        E = floyd_warshall_with_efficiency(G, weight_key="cost", eff_key="efficiency")
        assert np.all(D <= E)
        assert 2.0 not in E
        assert 2.0 not in D

    def test_return_all(self):
        G = nx.path_graph(10)
        add_edata(G)
        C, B, D = floyd_warshall_with_efficiency(
            G, weight_key="weight", eff_key="efficiency", return_all=True
        )
        assert np.all(B["weight"] >= C)
        assert np.all(B["efficiency"] <= 1)


class TestSymPyAllPairsShortestPath:
    def test_simple_graph(self):
        G, nodelist = simple_graph()
        C = sympy_floyd_warshall(
            G,
            f="weight / eff",
            accumulators={"weight": "sum", "eff": "product"},
            nonedge={"weight": np.inf, "eff": 0.0},
            nodelist=nodelist,
        )
        compare(G, C, nodelist)

    def test_simple_graph_new_eq(self):
        G, nodelist = simple_graph()
        for n1, n2, edata in G.edges(data=True):
            edata["eff"] = 1.0 / edata["eff"]
        C = sympy_floyd_warshall(
            G,
            f="weight * eff",
            accumulators={"weight": "sum", "eff": "product"},
            nodelist=nodelist,
        )
        for n1, n2, edata in G.edges(data=True):
            edata["eff"] = 1.0 / edata["eff"]
        compare(G, C, nodelist)

    @pytest.mark.parametrize("n", list(range(8)))
    @pytest.mark.parametrize("repeat", range(3))
    def test_complete(self, n, repeat):
        G, nodelist = complete_graph(n)
        C = sympy_floyd_warshall(
            G,
            f="weight / eff",
            accumulators={"weight": "sum", "eff": "product"},
            nonedge={"weight": np.inf, "eff": 0.0},
            nodelist=nodelist,
        )
        compare(G, C, nodelist)


class TestDijkstras:
    def test_simple_path(self):
        G = nx.path_graph(4)

        add_data(G, 0, 1, 100, 0.5)
        add_data(G, 1, 2, 100, 0.5)
        add_data(G, 2, 3, 100, 0.5)

        path_length, path = sympy_multisource_dijkstras(
            G, [0], "weight / eff", accumulators={"eff": "product"}
        )

        assert path_length[0] == 0.0
        assert path_length[1] == 200.0
        assert path_length[2] == 100 * 2 / 0.5 ** 2
        assert path_length[3] == 100 * 3 / 0.5 ** 3

    def test_simple_path(self):
        G = nx.path_graph(4)
        G.add_edge(0, 1, weight=100, eff=0.5)
        G.add_edge(1, 2, weight=100, eff=0.5)
        G.add_edge(2, 3, weight=100, eff=0.5)
        G.add_edge(3, 4, weight=100, eff=0.5)
        G.add_edge(0, 5, weight=200, eff=0.75)
        G.add_edge(5, 6, weight=200, eff=0.75)
        G.add_edge(6, 7, weight=200, eff=0.75)
        G.add_edge(7, 4, weight=200, eff=0.75)

        path_length, path = sympy_multisource_dijkstras(
            G, [0], "weight / eff", accumulators={"eff": "product"}
        )

        assert path_length[0] == 0.0
        assert path_length[1] == 100 * 1 / 0.5 ** 1
        assert path_length[2] == 100 * 2 / 0.5 ** 2
        assert path_length[3] == 100 * 3 / 0.5 ** 3
        assert path_length[4] == 200 * 4 / 0.75 ** 4
        assert path_length[5] == 200 * 1 / 0.75 ** 1
        assert path_length[6] == 200 * 2 / 0.75 ** 2
        assert path_length[7] == 200 * 3 / 0.75 ** 3

        assert path[3] == [0, 1, 2, 3]
        assert path[4] == [0, 5, 6, 7, 4]
        assert path[7] == [0, 5, 6, 7]

    def test_simple_path_source_to_target(self):
        G = nx.path_graph(4)
        G.add_edge(0, 1, weight=100, eff=0.5)
        G.add_edge(1, 2, weight=100, eff=0.5)
        G.add_edge(2, 3, weight=100, eff=0.5)
        G.add_edge(3, 4, weight=100, eff=0.5)
        G.add_edge(0, 5, weight=200, eff=0.75)
        G.add_edge(5, 6, weight=200, eff=0.75)
        G.add_edge(6, 7, weight=200, eff=0.75)
        G.add_edge(7, 4, weight=200, eff=0.75)

        path_length, path = sympy_dijkstras(
            g=G, source=0, target=4, f="weight / eff", accumulators={"eff": "product"}
        )
        assert path_length == 200 * 4 / 0.75 ** 4
        assert path == [0, 5, 6, 7, 4]


class TestBenchmarks:
    def test_dijkstras(self, benchmark):
        G, nodelist = complete_graph(200, create_using=nx.DiGraph)
        benchmark(
            floyd_warshall_with_efficiency,
            G,
            nodelist=nodelist,
            weight_key="weight",
            eff_key="efficiency",
        )
