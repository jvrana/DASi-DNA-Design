import networkx as nx
import numpy as np
from .utils import divide, select_from_arrs, replace_nan_with_inf
from sympy import sympify, lambdify
from collections import OrderedDict
from .exceptions import TerrariumNetworkxError

# TODO: cycle finder using *almost cycles*
# TODO: test directed vs undirected
# TODO: test multigraph


# def _archived_floyd_warhsall_with_efficiency(
#     G, weight="weight", efficiency="eff", nodelist=None, return_all=False, dtype=None
# ):
#     """
#     Computes the shortest path between all pairs using the cost function: SUM(w) / PROD(e)
#
#     Warning: This is guaranteed to return precise values due to floating point rounding errors.
#
#     :param G:
#     :param weight:
#     :param efficiency:
#     :param nodelist:
#     :return:
#     """
#     if dtype is None:
#         dtype = np.float64
#     A = nx.to_numpy_matrix(
#         G,
#         nodelist=nodelist,
#         multigraph_weight=min,
#         weight=weight,
#         nonedge=np.inf,
#         dtype=dtype,
#     )
#     B = nx.to_numpy_matrix(
#         G,
#         nodelist=nodelist,
#         multigraph_weight=min,
#         weight=efficiency,
#         nonedge=0.0,
#         dtype=dtype,
#     )
#
#     n, m = A.shape
#     I = np.identity(n)
#     A[I == 1] = 0  # diagonal elements should be zero
#     B[I == 1] = 1
#
#     for i in np.arange(n):
#         # get weight and efficiency of path using node 'i'
#         A_part = A[i, :] + A[:, i]
#         B_part = np.multiply(B[i, :], B[:, i])
#
#         # get total cost
#         C = divide(A, B)
#         C_part = divide(A_part, B_part)
#
#         # update
#         A = np.asmatrix(select_from_arrs(A, A_part, C < C_part))
#         B = np.asmatrix(select_from_arrs(B, B_part, C < C_part))
#
#     C = divide(A, B)
#     if return_all:
#         return C, A, B
#     else:
#         return C


PRODUCT = "product"
SUM = "sum"


def sympy_floyd_warshall(
    G,
    f,
    accumulators: dict,
    nonedge=None,
    nodelist=None,
    multigraph_weight=None,
    identity_subs=None,
    return_all=False,
    dtype=None,
):
    """
    Implementation of all pairs shortest path length using a modified floyd-warshall algorithm with arbitrary path
    length functions. The following path length function is valid:

    $$
    C = \frac{\sum_{i}^{n}{a_i}}{\prod_{i}^{n}{b_i}}
    $$

    Where $a_i$ and $b_i$ is the weight 'a' and 'b' of the *ith* edge in the path respectively.
    $\sum_{i}^{n}{a_i}$ is the accumulated sum of weights 'a' and $\prod_{i}^{n}{b_i}$
     is the accumulated product of weights 'b'. Arbitrarily complex path functions with arbitrary numbers of weights
     ($a, b, c,...$) can be used in the algorithm.

    Because arbitrary functions are used, the shortest path between ij and jk does not necessarily mean the shortest
    path nodes ijk is the concatenation of these two paths. In other words, for the shortest path $p_{ik}$ between
    nodes $i$ $j$ and $k$:

    $$
    p_{ij} + p_{jk} \neq p_{ijk}
    $$

    This means a predecessor dictionary cannot be used to reconstruct the paths easily. To do so, use the modified
    Dijkstra's algorithm below which handles arbitrary path functions.

    :param accumulators: if 'sum' or missing, apply $\sum$, if 'product' apply $\prod$ to edge weights
    :param nonedge:
    :param nodelist:
    :param multigraph_weight:
    :param identity_subs:
    :param return_all:
    :param dtype:
    :return:
    """
    if dtype is None:
        dtype = np.float64

    if identity_subs is None:
        identity_subs = {}

    if multigraph_weight is None:
        multigraph_weight = {}

    if nonedge is None:
        nonedge = {}

    expr = sympify(f)
    symbols = tuple(expr.free_symbols)
    func = lambdify(symbols, expr)

    matrix_dict = OrderedDict()

    for sym in symbols:
        matrix_dict[sym.name] = nx.to_numpy_matrix(
            G,
            nodelist=nodelist,
            multigraph_weight=multigraph_weight.get(sym.name, min),
            weight=sym.name,
            nonedge=nonedge.get(sym.name, np.inf),
            dtype=dtype,
        )

    n, m = list(matrix_dict.values())[0].shape

    # replace diagonals
    I = np.identity(n)
    for key, matrix in matrix_dict.items():
        # set accumulators
        if accumulators.get(key, SUM) == SUM:
            d = 0.0
        elif accumulators[key] == PRODUCT:
            d = 1.0
        else:
            raise TerrariumNetworkxError(
                "Accumulator key {} must either by '{}' or '{}' or a callable with two "
                "arguments ('M' a numpy matrix and 'i' a node index as an int)".format(
                    key
                )
            )

        # set diagonal
        matrix[I == 1] = identity_subs.get(key, d)

    for i in np.arange(n):
        # get costs if using node 'i'
        parts_dict = OrderedDict()
        for key, M in matrix_dict.items():
            M = matrix_dict[key]
            if accumulators.get(key, SUM) == SUM:
                parts_dict[key] = M[i, :] + M[:, i]
            elif accumulators[key] == PRODUCT:
                parts_dict[key] = np.multiply(M[i, :], M[:, i])
            else:
                raise TerrariumNetworkxError(
                    "Key '{}' not in accumulator dictionary. Options are '{}' or '{}'".format(
                        PRODUCT, SUM
                    )
                )

        # get total cost
        m_arr = [np.asarray(m) for m in matrix_dict.values()]
        p_arr = [np.asarray(m) for m in parts_dict.values()]
        C = replace_nan_with_inf(func(*m_arr))
        C_part = replace_nan_with_inf(func(*p_arr))

        # update
        for key, M in matrix_dict.items():
            part = parts_dict[key]
            c = C > C_part

            # TODO: is this if np.any(c) correct?
            if np.any(c):
                matrix_dict[key] = np.asmatrix(np.choose(c, (M, part)))

    m_arr = [np.asarray(m) for m in matrix_dict.values()]
    C = replace_nan_with_inf(func(*m_arr))
    if return_all:
        return C, matrix_dict
    return C


def floyd_warshall_with_efficiency(
    G, weight, efficiency, nodelist=None, return_all=False, dtype=None
):
    """
    Computes the shortest path between all pairs using the cost function: SUM(w) / PROD(e)

    Warning: This is guaranteed to return precise values due to floating point rounding errors.

    :param G:
    :param weight:
    :param efficiency:
    :param nodelist:
    :return:
    """
    f = "{} / {}".format(weight, efficiency)
    return sympy_floyd_warshall(
        G,
        f,
        accumulators={weight: SUM, efficiency: PRODUCT},
        nonedge={weight: np.inf, efficiency: 0.0},
        nodelist=nodelist,
        return_all=return_all,
        dtype=dtype,
    )
