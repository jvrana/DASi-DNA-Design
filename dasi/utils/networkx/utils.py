import numpy as np
from itertools import product, zip_longest
import networkx as nx


def select_from_arrs(A, B, condition):
    """Returns an ndarray of same shape as A and B, selecting elements from
    either A or B according to the condition. Somehow, this is 2X faster than
    whatever is implemented in np.choose."""
    a = np.asarray(A).ravel()
    b = np.asarray(B).ravel()
    c = np.asarray(condition).ravel()
    d = c * len(a) + np.arange(len(b))
    e = np.hstack((a, b))
    return e[d].reshape(A.shape)


def replace_nan_with_inf(m):
    m[np.isnan(m)] = np.inf
    return m


def divide(mata, matb):
    """Divides two matricies, replacing NaN with np.inf"""
    matc = mata / matb
    replace_nan_with_inf(matc)
    return matc


def find_all_min_paths(G, lim=None):
    """Explicitly finds minimum paths according to the cost function sum(weight) / prod(efficiency)"""
    min_path_dict = {}
    n_paths = 0
    for n1, n2 in product(G.nodes(), repeat=2):
        if lim and n_paths >= lim:
            break
        min_path = []
        min_cost = np.inf
        paths = nx.all_simple_paths(G, source=n1, target=n2)

        for path in map(nx.utils.pairwise, paths):
            path = list(path)
            weight = 0.0
            eff = 1.0
            for pair in path:
                edata = G.get_edge_data(*pair)
                weight += edata["weight"]
                eff *= edata["eff"]
            cost = weight / eff
            if cost < min_cost:
                min_path = path
                min_cost = cost
        min_path_dict.setdefault(n1, {})
        if n1 == n2:
            min_cost = 0.0
            min_path = []
        min_path_dict[n1][n2] = {"cost": min_cost, "path": min_path}
        n_paths += 1
    return min_path_dict


def compare_matrix_to_true_paths(mat, path_dict, nodelist):
    node_to_index = {n: i for i, n in enumerate(nodelist)}
    for n1 in path_dict:
        for n2, edata in path_dict[n1].items():
            u = node_to_index[n1]
            v = node_to_index[n2]

            min_cost = edata["cost"]
            est_cost = mat[u, v]
            if min_cost != est_cost:
                status = "ERROR {}".format(min_cost - est_cost)
            else:
                status = ""
            print(
                "{} -> {:5}: {:10.2f} {:13.2f} {:25}".format(
                    u, v, min_cost, est_cost, status
                )
            )


def min_edit_distance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for index2, char2 in enumerate(s2):
        newDistances = [index2 + 1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(
                    1
                    + min((distances[index1], distances[index1 + 1], newDistances[-1]))
                )
        distances = newDistances
    return distances[-1]


def compare_path_rank_order(mat, path_dict, nodelist, verbose=False):
    # get sorted index pairs representing end to end shortest paths
    indices = np.asarray(mat).flatten().argsort().flatten()
    index_pairs = np.dstack(np.unravel_index(indices, mat.shape))
    y = index_pairs.reshape(index_pairs.shape[1], index_pairs.shape[2])
    y = [(mat[u, v], u, v) for u, v in y]
    matrix_shortest_index_pairs = [(u, v) for c, u, v in sorted(y)]

    x = []
    node_to_index = {n: i for i, n in enumerate(nodelist)}
    for n1 in path_dict:
        for n2, edata in path_dict[n1].items():
            u = node_to_index[n1]
            v = node_to_index[n2]
            x.append((edata["cost"], u, v))
    true_shortest_index_pairs = [(u, v) for c, u, v in sorted(x)]

    edit_dist = min_edit_distance(
        true_shortest_index_pairs, matrix_shortest_index_pairs
    )
    print("Distance: {}".format(edit_dist))
    if verbose:
        for true_pair, mat_pair in zip_longest(
            true_shortest_index_pairs, matrix_shortest_index_pairs
        ):
            print("{}  {}".format(true_pair, mat_pair))


# TODO: move to networkx utils
def sort_cycle(arr, key=None):
    """Sort a cyclic array, maintaining order"""
    if key is None:
        arr_with_i = sorted([(x, i) for i, x in enumerate(arr)])
    else:
        arr_with_i = sorted([(key(x), i) for i, x in enumerate(arr)])
    i = arr_with_i[0][1]
    return arr[i:] + arr[:i]