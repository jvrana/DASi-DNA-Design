from typing import Dict
from typing import Hashable
from typing import List
from typing import Tuple

import networkx as nx
import numpy as np

from dasi.exceptions import DasiDesignException
from dasi.exceptions import DASiWarning
from dasi.models import AssemblyNode
from dasi.utils import sort_with_keys
from dasi.utils.networkx import sympy_floyd_warshall
from dasi.utils.networkx import sympy_multipoint_shortest_path
from warnings import warn

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


def only_long(w: np.ndarray, query_length: int, nodelist: List[Hashable]):
    # filter array based on nodes that span longer than the query
    index = np.array([n.index for n in nodelist]).reshape(1, -1)
    length = index - index.T
    w[np.where(length < query_length)] = np.inf
    return w


def only_short(w: np.ndarray, query_length: int, nodelist: List[Hashable]):
    # filter array based on nodes that span longer than the query
    index = np.array([n.index for n in nodelist]).reshape(1, -1)
    length = index - index.T
    w[np.where(length >= query_length)] = np.inf
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
        try:
            path = _multinode_to_shortest_path(graph, nodes, cyclic)
            if path not in unique_cyclic_paths:
                unique_cyclic_paths.append(path)
        except nx.NetworkXNoPath:
            pass
    return unique_cyclic_paths


def _get_closing_edge_indices(
    nodelist: List[AssemblyNode], query_length: int, matrix_dict
) -> Tuple[np.ndarray, np.ndarray]:
    node_to_i = {n: i for i, n in enumerate(nodelist)}
    src, dest = [], []

    # closing edge indices
    for n in nodelist:
        i1 = node_to_i[n]
        n2 = AssemblyNode(n.index + query_length, *list(n)[1:])
        if n2 in node_to_i:
            i2 = node_to_i[n2]
            v = matrix_dict["material"][i1, i2]
            if not np.isinf(v):
                src.append(i1)
                dest.append(i2)
    index = (np.array(src), np.array(dest))
    return index


def _get_closure_matrix(
    mat: np.ndarray,
    eff: np.ndarray,
    nodelist: List[AssemblyNode],
    query_length: int,
    matrix_dict: dict,
) -> Tuple[np.ndarray, np.ndarray]:
    indices = _get_closing_edge_indices(nodelist, query_length, matrix_dict)
    m = np.full_like(mat, np.inf)
    e = np.full_like(eff, 0.0)
    if len(indices[0]):
        m[indices] = mat[indices]
        e[indices] = eff[indices]
    return m, e


def cyclic_matrix(
    matrix_dict: Dict[str, np.matrix], nodelist: List[AssemblyNode], query_length: int
) -> np.ndarray:
    m1 = np.array(matrix_dict["material"])
    e1 = np.array(matrix_dict["efficiency"])
    m2, e2 = _get_closure_matrix(m1, e1, nodelist, query_length, matrix_dict)

    m = m1 + m2.min(1)
    e = e1 * e2.max(1)

    w = np.divide(m, e)
    fill_diag_inf(w)
    return w


def optimize_graph(
    graph: nx.DiGraph, query_length: int, cyclic: bool, n_paths: int
) -> Tuple[List[List[tuple]], List[float]]:
    # get ordered nodelist and nodekeys
    nodelist, nodekeys = sort_with_keys(list(graph.nodes()), key=lambda x: x[0])

    # 2D matrix of 'efficiency' and 'material' costs
    weight_matrix, matrix_dict, ori_matrix_dict = sympy_floyd_warshall(
        graph,
        f=path_length_config["f"],
        accumulators=path_length_config["accumulators"],
        nodelist=nodelist,
        dtype=np.float64,
        return_all=True,
    )
    matrix_dict = {k: np.array(v) for k, v in matrix_dict.items()}

    if cyclic:
        weight_matrix = cyclic_matrix(matrix_dict, nodelist, query_length)
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
    paths = _nodes_to_fullpaths(graph, trimmed_nodes, cyclic, n_paths=n_paths)

    if len(paths) < n_paths:
        warn("Number of paths returned is less than requested paths {} < {}".format(
                len(paths), n_paths), DASiWarning)

    _check_paths(paths)
    return paths, trimmed_costs
