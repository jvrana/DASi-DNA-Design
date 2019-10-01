"""Utilities

.. module:: dasi.utils

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    networkx
    npdf
    async_wrapper
    region
    span
"""

import bisect
from .region import Region
import functools
from .npdf import NumpyDataFrame, NumpyDataFrameException
from .networkx.utils import sort_cycle
from .networkx.shortest_path import multipoint_shortest_path


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
    if data["strand"] == 1 and data["start"] == 1 and data["raw_end"] == data["length"]:
        return True
    elif (
        data["strand"] == -1
        and data["raw_end"] == 1
        and data["start"] == data["length"]
    ):
        return True


def partialclass(cls, *args, **kwds):
    class PartialClass(cls):
        __init__ = functools.partialmethod(cls.__init__, *args, **kwds)

    return PartialClass
