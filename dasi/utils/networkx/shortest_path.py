from heapq import heappop
from heapq import heappush
from itertools import count
from typing import List

import networkx as nx
import numpy as np
from more_itertools import pairwise
from sympy import lambdify
from sympy import sympify

from .exceptions import NetworkxUtilsException
from .utils import sort_cycle


def _weight_function(v, u, e, k):
    if e is None:
        return None
    return e[k]


def sympy_multisource_dijkstras(
    g, sources, f, accumulators=None, init=None, target=None, cutoff=None
):
    if sources is None or not len(sources):
        raise ValueError("sources must not be empty")
    if target in sources:
        return (0, [target])
    paths = {source: [source] for source in sources}  # dictionary of paths
    dist = _multisource_dijkstra(
        g,
        sources,
        f,
        target=target,
        accumulators=accumulators,
        init=init,
        paths=paths,
        cutoff=cutoff,
    )
    if target is None:
        return (dist, paths)
    try:
        return (dist[target], paths[target])
    except KeyError:
        raise nx.NetworkXNoPath("No path to {}.".format(target))


def sympy_dijkstras(
    g, source, f, target=None, accumulators=None, init=None, cutoff=None
):
    """Computes the shortest path distance and path for a graph using an
    arbitrary function.

    :param g:
    :param source:
    :param f:
    :param target:
    :param accumulators:
    :param init:
    :param cutoff:
    :return:
    """
    dist, path = sympy_multisource_dijkstras(
        g,
        [source],
        f,
        target=target,
        accumulators=accumulators,
        init=init,
        cutoff=cutoff,
    )
    return dist, path


def _multisource_dijkstra(
    g,
    sources,
    f,
    target=None,
    accumulators=None,
    init=None,
    cutoff=None,
    paths=None,
    pred=None,
):
    accumulators = accumulators or {}
    init = init or {}

    # successor dictionary
    g_succ = g._succ if g.is_directed() else g._adj

    # push/pop methods to use
    push = heappush
    pop = heappop

    # node to shortest distance, breakdown
    dist_parts = {}

    # node to shortest distance
    dist = {}

    # node to shortest distance seen
    seen = {}
    c = count()
    fringe = []

    # lambda function
    func = sympify(f)
    symbols = list([s.name for s in func.free_symbols])
    func = lambdify(symbols, func)

    # initial/default values for each symbol
    for sym in symbols:
        if accumulators.get(sym, "sum") == "sum":
            init.setdefault(sym, 0.0)
        elif accumulators.get(sym, "product"):
            init.setdefault(sym, 1.0)
        else:
            raise NetworkxUtilsException("Accumulator '{}' not recognized".format(sym))
    init = np.array([init[x] for x in symbols])

    # accumulator functions to each for each symbol
    # if 'product', return the accumulated product
    # else return the accumulated sum (default)
    accu_f = []
    for sym in symbols:
        if accumulators.get(sym, "sum") == "sum":
            accu_f.append(lambda x: np.sum(x))
        elif accumulators[sym] == "product":
            accu_f.append(lambda x: np.prod(x))

    # push the initial values for the sources
    for source in sources:
        if source not in g:
            raise nx.NodeNotFound("Source {} is not in G".format(source))
        seen[source] = func(*init)
        push(fringe, (0, next(c), source, init))

    # modified dijkstra's
    while fringe:
        (_, _, v, d) = pop(fringe)
        # d np.array of values
        if v in dist_parts:
            continue  # already searched this node
        dist_parts[v] = d
        dist[v] = func(*d)
        if v == target:
            break
        for u, e in g_succ[v].items():
            # vu cost breakdown for  each symbol
            costs = np.array([_weight_function(v, u, e, sym) for sym in symbols])

            # vu_dist break down using accumulating function
            x = np.stack([dist_parts[v], costs])
            vu_dist_parts = []
            for i, _x in enumerate(x.T):
                vu_dist_parts.append(accu_f[i](_x))
            vu_dist_parts = np.array(vu_dist_parts)
            vu_dist = func(*vu_dist_parts)

            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            if u in dist_parts:
                if vu_dist < dist[u]:
                    raise ValueError("Contradictory paths found:", "negative weights?")
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = func(*vu_dist_parts)
                push(fringe, (vu_dist, next(c), u, vu_dist_parts))
                if paths is not None:
                    paths[u] = paths[v] + [u]
                if pred is not None:
                    pred[u] = v
            elif vu_dist == seen[u]:
                if pred is not None:
                    pred[u].append(v)
                    fringe = []
    return dist


def multipoint_shortest_path(
    graph: nx.DiGraph,
    nodes: List[str],
    weight_key: str,
    cyclic=False,
    cyclic_sort_key=None,
):
    """Return shortest path through nodes. If cyclic, will return the cycle
    sorted with the 'lowest' node at index 0. Self cycles are not supported.

    :param graph: the graph
    :param nodes: list of nodes to find path
    :param weight_key: weight key
    :param cyclic: whether the path is cyclic
    :param cyclic_sort_key: the key function to use to sort the cycle (if cyclic)
    :return:
    """
    if cyclic_sort_key and not cyclic:
        raise ValueError("cyclic_sort_key was provided but 'cyclic' was False.")
    full_path = []
    if cyclic:
        nodes = nodes + nodes[:1]
    for n1, n2 in pairwise(nodes):
        path = nx.shortest_path(graph, n1, n2, weight=weight_key)
        full_path += path[:-1]
    if not cyclic:
        full_path.append(nodes[-1])
    if cyclic:
        return sort_cycle(full_path, cyclic_sort_key)
    else:
        return full_path


def sympy_multipoint_shortest_path(
    graph: nx.DiGraph,
    nodes: List[str],
    f: str,
    accumulators: dict,
    init=None,
    cutoff=None,
    cyclic=False,
    cyclic_sort_key=None,
):
    if cyclic_sort_key and not cyclic:
        raise ValueError("cyclic_sort_key was provided but 'cyclic' was False.")
    full_path = []
    full_path_length = 0.0
    if cyclic:
        nodes = nodes + nodes[:1]
    for n1, n2 in pairwise(nodes):
        path_length, path = sympy_dijkstras(
            graph,
            f=f,
            source=n1,
            target=n2,
            accumulators=accumulators,
            init=init,
            cutoff=cutoff,
        )
        full_path_length += path_length
        full_path += path[:-1]
    if not cyclic:
        full_path.append(nodes[-1])
    if cyclic:
        return full_path_length, sort_cycle(full_path, cyclic_sort_key)
    else:
        return full_path_length, full_path
