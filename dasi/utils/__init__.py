"""Utilities

.. module:: design

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    async_wrapper
    region
    span
"""

import bisect
from .region import Region
from .span import Span
import functools
from .async_wrapper import make_async
from more_itertools import pairwise
import networkx as nx
from typing import List


def sort_with_keys(a, key):
    s = sorted(a, key=key)
    keys = [key(x) for x in s]
    return s, keys


def bisect_between(a, low, high):
    """
    Returns the start (inclusive) and end (exclusive) indices
    for a sorted array using a key function.

    :param a: sorted array (does not check)
    :param low: low key
    :param high: high key
    :param key: key function
    :return: tuple of start (inclusive) and end (exclusive) indices
    """
    i = bisect.bisect_left(a, low)
    j = bisect.bisect_right(sorted(a[i:]), high)
    return i, j + i


def bisect_slice_between(a, keys, low, high):
    i, j = bisect_between(keys, low, high)
    return a[i:j]


def perfect_subject(data):
    """Determine whether a blast result consumes the entire subject."""
    if data["strand"] == 1 and data["start"] == 1 and data["end"] == data["length"]:
        return True
    elif data["strand"] == -1 and data["end"] == 1 and data["start"] == data["length"]:
        return True


def partialclass(cls, *args, **kwds):
    class PartialClass(cls):
        __init__ = functools.partialmethod(cls.__init__, *args, **kwds)

    return PartialClass


def sort_cycle(arr, key=None):
    """Sort a cyclic array, maintaining order"""
    if key is None:
        arr_with_i = sorted([(x, i) for i, x in enumerate(arr)])
    else:
        arr_with_i = sorted([(key(x), i) for i, x in enumerate(arr)])
    i = arr_with_i[0][1]
    return arr[i:] + arr[:i]


def multipoint_shortest_path(
    graph: nx.DiGraph,
    nodes: List[str],
    weight_key: str,
    cyclic=False,
    cyclic_sort_key=None,
):
    """
    Return shortest path through nodes. If cyclic, will return the cycle sorted with the
    'lowest' node at index 0. Self cycles are not supported

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
