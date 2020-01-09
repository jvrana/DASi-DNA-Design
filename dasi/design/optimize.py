from collections import OrderedDict
from typing import List
from typing import Tuple

import networkx as nx
import numpy as np

from dasi.exceptions import DasiDesignException
from dasi.exceptions import DASiWarning
from dasi.utils import sort_with_keys
from dasi.utils.networkx import sympy_floyd_warshall
from dasi.utils.networkx import sympy_multipoint_shortest_path
from dasi.utils.networkx.algorithms import accumulate_helper
from dasi.utils.networkx.algorithms import str_to_symbols_and_func
from dasi.utils.networkx.utils import replace_nan_with_inf

# definition of how to compute path length
# c = SUM(m) / PRODUCT(e), where m and e are arrays of attributes 'material'
# and 'efficiency' for a given path

path_length_config = {
    "f": "material / efficiency ",
    "accumulators": {"efficiency": "product"},
}


def _check_paths(paths: List[List[tuple]]) -> None:
    """Validates a path to check for duplicate nodes.

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


def argmin(w):
    index = np.argsort(w.ravel())
    return np.unravel_index(index, w.shape)


def fill_diag(w, fillval):
    iden = np.identity(len(w))
    w[np.where(iden == 1)] = fillval
    return w


def fill_diag_inf(w):
    return fill_diag(w, np.inf)


def index_slice(indices, arr):
    """From tuple of indices, return zipped items."""
    unzipped = []
    for i in indices:
        unzipped.append([arr[_i] for _i in i])
    return list(zip(*unzipped))


# TODO: is this necessary?
def only_ab(w, nodelist):
    b_array = np.array([n.type == "B" for n in nodelist]).reshape(1, -1)
    w[np.where(~np.logical_xor(b_array, b_array.T))] = np.inf
    return w


def only_long(w, query_length, nodelist):
    index = np.array([n.index for n in nodelist]).reshape(1, -1)
    length = index - index.T
    w[np.where(length < query_length)] = np.inf
    return w


def _multinode_to_shortest_path(
    graph, nodes, cyclic, return_length=False
) -> List[tuple]:
    """Estimate the shortest path that touches the specified nodes."""
    path_length, path = sympy_multipoint_shortest_path(
        graph,
        nodes,
        f=path_length_config["f"],
        accumulators=path_length_config["accumulators"],
        cyclic=cyclic,
    )
    if return_length:
        return path, path_length
    return path


def _nodes_to_fullpaths(
    graph: nx.DiGraph, cycle_endpoints: Tuple, cyclic: bool, n_paths=None
) -> List[List[Tuple]]:
    """Recover full paths from cycle endpoints.

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


def optimize_graph(
    graph: nx.DiGraph, query_length: int, cyclic: bool, n_paths: int
) -> Tuple[List[List[tuple]], List[float]]:
    nodelist, nodekeys = sort_with_keys(list(graph.nodes()), key=lambda x: x[0])
    weight_matrix, matrix_dict, ori_matrix_dict = sympy_floyd_warshall(
        graph,
        f=path_length_config["f"],
        accumulators=path_length_config["accumulators"],
        nodelist=nodelist,
        dtype=np.float64,
        return_all=True,
    )

    if cyclic:
        # add the closing edge
        closed_matrix_dict = OrderedDict({k: v.copy() for k, v in matrix_dict.items()})

        # # fold at length
        # fold_at_length = []
        # for n in nodelist:
        #     if n.index > query_length:
        #         n = AssemblyNode(n.index - query_length, *list(n)[1:])
        #     fold_at_length.append(node_to_i[n])

        for k, v in closed_matrix_dict.items():
            m1 = matrix_dict[k].copy()
            m1 = only_long(m1, query_length, nodelist)

            # TODO: do we need to fold this?
            m2 = ori_matrix_dict[k][:, :].T
            # m2 = ori_matrix_dict[k][:, fold_at_length].T
            closed_matrix_dict[k] = accumulate_helper(
                path_length_config["accumulators"].get(k, "sum"), m1, m2
            )

        symbols, func = str_to_symbols_and_func(path_length_config["f"])
        closed_wmatrix = func(*[np.asarray(m) for m in closed_matrix_dict.values()])
        closed_wmatrix = replace_nan_with_inf(closed_wmatrix)

        fill_diag_inf(closed_wmatrix)
        # only_ab(closed_wmatrix, nodelist)
        weight_matrix = closed_wmatrix
    else:
        raise NotImplementedError("Linear assemblies not yet implemented.")

    min_index = tuple([i[:n_paths] for i in argmin(weight_matrix)])
    costs = weight_matrix[min_index]
    costs = [c for c in costs if c != np.inf]
    a_nodes = [nodelist[i] for i in min_index[0]]
    b_nodes = [nodelist[i] for i in min_index[1]]
    nodes = list(zip(a_nodes, b_nodes))

    nodes_and_costs = list(zip(nodes, costs))

    if nodes_and_costs:
        trimmed_nodes, trimmed_costs = zip(*nodes_and_costs)
    else:
        trimmed_nodes, trimmed_costs = [], []
    paths = _nodes_to_fullpaths(graph, trimmed_nodes, False, n_paths=n_paths)

    if len(paths) < n_paths:
        DASiWarning(
            "Number of paths returned is less than requested paths {} < {}".format(
                len(paths), n_paths
            )
        )

    _check_paths(paths)
    return paths, trimmed_costs
