r"""
Utilities (:mod:`dasi.utils`)
=============================

This module provide various utility functions.

.. currentmodule:: dasi.utils

Utility modules
---------------

.. autosummary::
    :toctree: generated/

    npdf
    region
    sequence_design

Networkx utilities
------------------

.. currentmodule:: dasi.utils.networkx

Specialized networkx algorithms for path and cycle finding.

.. autosummary::
    :toctree: generated/

    algorithsm
    exceptions
    shortest_path
    utils

"""
import bisect
from typing import Any
from typing import Callable
from typing import Iterable
from typing import List
from typing import Tuple

from .networkx.shortest_path import multipoint_shortest_path
from .networkx.utils import sort_cycle
from .npdf import NumpyDataFrame
from .npdf import NumpyDataFrameException
from .region import Region


def sort_with_keys(a: Iterable[Any], key: Callable) -> Tuple[List, List]:
    """Sort an iterable, returning both the sorted array and the sorted keys.

    :param a: the iterable
    :param key: key function to use for sorting
    :return:
    """
    s = sorted(a, key=key)
    keys = [key(x) for x in s]
    return s, keys


def bisect_between(a: Iterable, low: Any, high: Any) -> Tuple[int, int]:
    """Returns the start (inclusive) and end (exclusive) indices for a sorted
    array using a key function.

    :param a: sorted array (does not check)
    :param low: low key
    :param high: high key
    :return: tuple of start (inclusive) and end (exclusive) indices
    """
    i = bisect.bisect_left(a, low)
    j = bisect.bisect_right(a[i:], high)
    return i, j + i


def bisect_slice_between(a: Iterable, keys: Iterable, low: Any, high: Any) -> Iterable:
    """Slice the iterable using inclusive bisection. Assumes both the iterable
    and keys are sorted. Bisect at specified `low` and `high`.

    :param a: pre-sorted iterable to slice
    :param keys: pre-sorted keys to bisect
    :param low: low key
    :param high: high key
    :return: sliced iterable
    """
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


def prep_df(df):
    colnames = [
        "DESIGN_ID",
        "DESIGN_KEY",
        "ASSEMBLY_ID",
        "REACTION_ID",
        "NAME",
        "TYPE",
        "KEY",
        "ROLE",
        "REGION",
        "SEQUENCE",
        "LENGTH",
        "META",
    ]
    df.columns = colnames
    df.sort_values(by=["TYPE", "DESIGN_ID", "REACTION_ID", "ASSEMBLY_ID", "ROLE"])
    return df


def group_by(arr: List[Any], key: Callable):
    """Group a list by some key."""
    grouped = {}
    for x in arr:
        k = key(x)
        grouped.setdefault(k, list())
        grouped[k].append(x)
    return grouped
