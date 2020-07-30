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
import inspect
from copy import deepcopy
from datetime import datetime
from functools import wraps
from itertools import tee
from typing import Any
from typing import Callable
from typing import Generator
from typing import Iterable
from typing import List
from typing import Tuple
from typing import TypeVar
from typing import Union

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


def now():
    return datetime.now()


def log_times(key: str = None, class_attribute: str = "_method_trace"):
    """wrapper for logging method run times for a class."""

    def wrapped(f):
        @wraps(f)
        def _wrapped(self, *args, **kwargs):

            if not hasattr(self, class_attribute):
                raise ValueError(
                    "Instance {} must have attribute '{}'".format(self, class_attribute)
                )
            elif not isinstance(getattr(self, class_attribute), dict):
                raise ValueError(
                    "Attribute {} of {} must be a {}".format(
                        class_attribute, self, dict
                    )
                )
            t1 = now()
            result = f(self, *args, **kwargs)
            t2 = now()

            if key is None:
                use_key = f.__name__
            else:
                use_key = key
            getattr(self, class_attribute)[use_key] = (t1, t2)
            return result

        return _wrapped

    return wrapped


def fmt_datetime(t):
    return str(t)


def log_metadata(
    key: str = None,
    class_attribute: str = "_method_trace",
    additional_metadata: dict = None,
):
    """wrapper for logging method run times for a class."""

    def wrapped(f):
        @wraps(f)
        def _wrapped(self, *args, **kwargs):

            if not hasattr(self, class_attribute):
                raise ValueError(
                    "Instance {} must have attribute '{}'".format(self, class_attribute)
                )
            elif not isinstance(getattr(self, class_attribute), dict):
                raise ValueError(
                    "Attribute {} of {} must be a {}".format(
                        class_attribute, self, dict
                    )
                )
            t1 = now()
            result = f(self, *args, **kwargs)
            t2 = now()

            argspec = inspect.getfullargspec(f)

            copied_args = deepcopy(args)
            copied_kwargs = deepcopy(kwargs)
            argdict = dict(zip(argspec.args[1:], copied_args))
            argdict.update(copied_kwargs)

            metadata = {
                "__name__": f.__name__,
                "__spec__": str(argspec),
                "start": fmt_datetime(t1),
                "end": fmt_datetime(t2),
                "args": argdict,
            }
            if additional_metadata:
                metadata.update(additional_metadata)

            if key is None:
                use_key = f.__name__
            else:
                use_key = key
            getattr(self, class_attribute)[use_key] = metadata
            return result

        return _wrapped

    return wrapped


T = TypeVar("T")


def argsorted(
    arr: Iterable[T], key: Callable, return_items: bool = False
) -> Union[List[Tuple[int, T]], List[T]]:
    s = sorted(enumerate(tee(arr)[0]), key=lambda x: key(x[1]))
    if return_items:
        return s
    else:
        return [_s[0] for _s in s]


def lexsorted(keys: Iterable, target: Iterable[T], key: Callable) -> List[T]:
    sorted_indices = argsorted(enumerate(tee(keys)[0]), key=key)
    return [target[i] for i in sorted_indices]


def chunkify(arr: Iterable[T], chunk_size: int) -> Generator[List[T], None, None]:
    new_list = []
    for x in tee(arr, 1)[0]:
        new_list.append(x)
        if len(new_list) == chunk_size:
            yield new_list
            new_list = []
    if new_list:
        yield new_list
