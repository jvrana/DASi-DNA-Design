from collections import OrderedDict
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
from sympy import lambdify
from sympy import sympify

from .exceptions import TerrariumNetworkxError
from .utils import replace_nan_with_inf
from .utils import select_from_arrs

PRODUCT = "product"
SUM = "sum"


def sympy_floyd_warshall(
    g: Union[nx.DiGraph, nx.Graph],
    f: str,
    accumulators: dict,
    nonedge=None,
    nodelist=None,
    multigraph_weight=None,
    identity_subs=None,
    return_all=False,
    dtype=None,
) -> Union[np.ndarray, Tuple[np.ndarray, dict]]:
    """
    Implementation of all pairs shortest path length using a modified floyd-warshall
    algorithm with arbitrary path
    length functions. The following path length function is valid:

    $$
    C = \frac{\\sum_{i}^{n}{a_i}}{\\prod_{i}^{n}{b_i}}
    $$

    Where $a_i$ and $b_i$ is the weight 'a' and 'b' of the *ith* edge in the path
    respectively.
    $\\sum_{i}^{n}{a_i}$ is the accumulated sum of weights 'a' and $\\prod_{i}^{n}{b_i}$
     is the accumulated product of weights 'b'. Arbitrarily complex path functions with
      arbitrary numbers of weights
     ($a, b, c,...$) can be used in the algorithm.

    Because arbitrary functions are used, the shortest path between ij and jk does not
    necessarily mean the shortest
    path nodes ijk is the concatenation of these two paths. In other words, for the
    shortest path $p_{ik}$ between
    nodes $i$ $j$ and $k$:

    $$
    p_{ij} + p_{jk} \neq p_{ijk}
    $$

    This means a predecessor dictionary cannot be used to reconstruct the paths easily.
    To do so, use the modified
    Dijkstra's algorithm below which handles arbitrary path functions.

    :param g: the graph
    :param f: the function string that represents SymPy function to compute the weights
    :param accumulators: diciontary of symbols to accumulator functions (choose from
                         ["PRODUCT", "SUM"] to use for accumulation of weights through
                         a path. If missing "SUM" is used.
    :param nonedge: dictionary of symbol to value to use for nonedges
                    (e.g. {'weight': np.inf})
    :param nodelist: optional nodelist to use
    :param multigraph_weight: optional (default: min) function to use for multigraphs
    :param identity_subs: the dictionary of values to set along the diagonal axis
    :param return_all: if True, return both the resulting weight matrix and the
            individual components broken down by symbol strings.
    :param dtype: the dtype of the resulting np.ndarrays used and returned
    :return: either just the weight matrix or, if return_all is True, the weight_matrix
                and the dictionary of the weight components.
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
            g,
            nodelist=nodelist,
            multigraph_weight=multigraph_weight.get(sym.name, min),
            weight=sym.name,
            nonedge=nonedge.get(sym.name, np.inf),
            dtype=dtype,
        )

    n, m = list(matrix_dict.values())[0].shape

    # replace diagonals
    identity = np.identity(n)
    for key, matrix in matrix_dict.items():
        # set accumulators
        if accumulators.get(key, SUM) == SUM:
            d = 0.0
        elif accumulators[key] == PRODUCT:
            d = 1.0
        else:
            raise TerrariumNetworkxError(
                "Accumulator key {} must either be '{}' or '{}' or a callable with two "
                "arguments ('M' a numpy matrix and 'i' a node index as an int)".format(
                    key, SUM, PRODUCT
                )
            )

        # set diagonal
        matrix[identity == 1] = identity_subs.get(key, d)

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
                        key, PRODUCT, SUM
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

            if np.any(c):
                assert M.shape == part.shape
                # matrix_dict[key] = np.asmatrix(np.choose(c, (M, part)))
                matrix_dict[key] = np.asmatrix(select_from_arrs(part, M, C < C_part))

    m_arr = [np.asarray(m) for m in matrix_dict.values()]
    C = replace_nan_with_inf(func(*m_arr))
    if return_all:
        return C, matrix_dict
    return C


def floyd_warshall_with_efficiency(
    g: Union[nx.DiGraph, nx.Graph],
    weight_key: str,
    eff_key: str,
    nodelist=None,
    return_all=False,
    dtype=None,
) -> np.ndarray:
    """Computes the shortest path between all pairs using the cost function:
    SUM(w) / PROD(e)

    Warning: This is *not guaranteed* to return precise values due to floating point
    rounding errors.

    :param g: the graph
    :param weight_key: the weight key
    :param eff_key: the efficiency key
    :param nodelist: optional list of nodes
    :param return_all: whether to return the weight matrix and weight dictionary
    :param dtype: dtype of the np.array matrix
    :return: the weight matrix
    """
    f = "{} / {}".format(weight_key, eff_key)
    return sympy_floyd_warshall(
        g,
        f,
        accumulators={weight_key: SUM, eff_key: PRODUCT},
        nonedge={weight_key: np.inf, eff_key: 0.0},
        nodelist=nodelist,
        return_all=return_all,
        dtype=dtype,
    )
