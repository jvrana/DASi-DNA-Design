import bisect
from typing import List, Tuple

import networkx as nx
import numpy as np

from dasi.design.graph_builder import AssemblyGraphBuilder
from dasi.exceptions import DasiDesignException
from dasi.utils import sort_with_keys
from dasi.utils.networkx import sympy_floyd_warshall, sympy_multipoint_shortest_path
from multiprocessing import Pool

# definition of how to compute path length
# c = SUM(m) / PRODUCT(e), where m and e are arrays of attributes 'material' and 'efficiency' for a given path
path_length_config = {
    "f": "material / efficiency",
    "accumulators": {"efficiency": "product"},
}


def _check_paths(paths):
    """
    Validates a path to check for duplicate nodes.

    :param paths:
    :return:
    """
    invalid_paths = []
    for path in paths:
        lastseen = path[0][2]
        for p in path[1:]:
            if p[2] == lastseen:
                invalid_paths.append(path)
                break
            lastseen = p[2]
    if invalid_paths:
        raise DasiDesignException(
            "There are {} invalid paths:\n{}\n...{} more".format(
                len(invalid_paths),
                "\n".join([str(x) for x in invalid_paths[:5]]),
                max(len(invalid_paths) - 5, 0),
            )
        )


def _all_pairs_shortest_path(graph, nodelist):
    return sympy_floyd_warshall(
        graph,
        f=path_length_config["f"],
        accumulators=path_length_config["accumulators"],
        nodelist=nodelist,
        dtype=np.float64,
    )


def _multinode_to_shortest_path(graph, nodes, cyclic):
    """Estimate the shortest path that touches the specified nodes."""
    path_length, path = sympy_multipoint_shortest_path(
        graph,
        nodes,
        f=path_length_config["f"],
        accumulators=path_length_config["accumulators"],
        cyclic=cyclic,
    )
    return path


def _collect_cycle_endpoints(graph: nx.DiGraph, length: int):
    """
    Use the floyd-warshall algorithm to compute the shortest path endpoints.

    :param graph: the networkx graph
    :param length: the size of the query sequence
    :return:
    """
    nodelist, nkeys = sort_with_keys(list(graph.nodes()), key=lambda x: x[0])
    node_to_i = {v: i for i, v in enumerate(nodelist)}
    weight_matrix = _all_pairs_shortest_path(graph, nodelist)
    endpoints = []

    def bisect_iterator(nodelist, nkeys):
        _i = bisect.bisect_right(nkeys, length)
        for i, A in enumerate(nodelist[:_i]):
            _j = bisect.bisect_left(nkeys, A.index + length)
            for B in nodelist[_j:]:
                if B.type == "B":
                    j = node_to_i[B]
                    yield i, j, A, B

    pair_iterator = bisect_iterator(nodelist, nkeys)
    for i, j, A, B in pair_iterator:

        # TODO: must include final edge
        a = weight_matrix[i, j]
        b = weight_matrix[j, i]
        if a != np.inf and b != np.inf:
            x = ((A, B), (a, b), a + b)
            endpoints.append(x)

    endpoints = sorted(endpoints, key=lambda x: (x[-1], x[0]))
    return endpoints


def _nodes_to_fullpaths(
    graph: nx.DiGraph, cycle_endpoints: Tuple, cyclic: bool, n_paths=None
) -> List[List[Tuple]]:
    """
    Recover full paths from cycle endpoints.

    :param graph:
    :param cycle_endpoints:
    :param n_paths:
    :return:
    """
    unique_cyclic_paths = []
    for nodes in cycle_endpoints:
        if n_paths is not None and len(unique_cyclic_paths) >= n_paths:
            break
        path = _multinode_to_shortest_path(graph, nodes, cyclic)
        if path not in unique_cyclic_paths:
            unique_cyclic_paths.append(path)
    return unique_cyclic_paths


def _collect_optimized_paths(
    graph: nx.DiGraph, length: int, cyclic: bool, n_paths=None
):
    """
    Collect minimum cycles or linear paths from a graph.

    :param graph: the networkx graph representing the assembly graph
    :param length: length of the query
    :param cyclic: whether to search for a cyclic assembly
    :param n_paths: maximum number of paths to return
    :return:
    """
    if cyclic:
        nodes = _collect_cycle_endpoints(graph, length=length)
    else:
        raise NotImplementedError("Linear assemblies are not yet implemented.")
    paths = _nodes_to_fullpaths(
        graph, [n[0] for n in nodes], cyclic=cyclic, n_paths=n_paths
    )
    _check_paths(paths)
    return paths


def optimize_graph(graph, query_length, cyclic, n_paths):
    """Optimize the graph associated with the specified query_key"""
    # self.logger.info("Optimizing {}".format(query_key))
    paths = _collect_optimized_paths(graph, query_length, cyclic, n_paths=n_paths)
    return paths


def _multiprocessing_optimize_graph(args):
    return optimize_graph(args[0], args[1], args[2], args[3])


def multiprocessing_optimize_graph(graphs, query_lengths, cyclics, n_paths, n_jobs):
    """Optimize graphs using multiprocessing"""
    args = [(g, q, c, n_paths) for g, q, c in zip(graphs, query_lengths, cyclics)]

    with Pool(processes=min(n_jobs, len(graphs))) as pool:  # start 4 worker processes
        paths = pool.map(_multiprocessing_optimize_graph, args)
    return paths


def assemble_graph(container, span_cost):
    """Build an assembly graph for a specified query."""
    container.expand(expand_overlaps=True, expand_primers=True)
    container.freeze()
    # group by query_regions
    # groups = container.groups()
    # self.logger.info(
    #     "Number of types: {}".format(len(container.groups_by_type))
    # )
    # self.logger.info("Number of groups: {}".format(len(groups)))
    # build assembly graph
    graph_builder = AssemblyGraphBuilder(container, span_cost=span_cost)
    G = graph_builder.build_assembly_graph()
    # self.logger.info("=== Assembly Graph ===")
    # self.logger.info(nx.info(G))
    assert G.number_of_edges()
    return G


def _multiprocessing_assemble_graph(arg):
    return assemble_graph(arg[0], arg[1])


def multiprocessing_assemble_graph(containers, span_cost, n_jobs):
    """Assemble graphs using multiprocessing"""
    args = [(container, span_cost) for container in containers]
    with Pool(
        processes=min(n_jobs, len(containers))
    ) as pool:  # start 4 worker processes
        graphs = pool.map(_multiprocessing_assemble_graph, args)
    return graphs
