import functools
from functools import lru_cache
from itertools import repeat
from itertools import zip_longest
from typing import Iterable
from typing import List
from typing import Optional
from typing import Union

import numpy as np

from dasi.utils.sequence.sequence_complexity import DNAStats


@lru_cache(254)
def cached_stats_cost(stats: DNAStats, i: int, j: int):
    return stats.cost(i, j)


def stats_cost_arr(stats: DNAStats, i: Iterable[int], j: Iterable[int]):
    for _i, _j in zip_longest(i, j, fillvalue="STOP"):
        if i == "STOP":
            raise ValueError("Iterables are different sizes.")
        elif j == "STOP":
            raise ValueError("Iterables are different sizes")
        yield cached_stats_cost(stats, _i, _j)


def get_partitions(stats: DNAStats, x: List[int], delta: int = 0):
    """Get partition costs for a list of partition indices."""
    x1 = np.array(x)
    if delta is None:
        x2 = x1
    else:
        x2 = x1 + delta
    c1 = np.array(list(stats_cost_arr(stats, repeat(None, len(x1)), x1)))
    c2 = np.array(list(stats_cost_arr(stats, x2, repeat(None, (len(x2))))))
    return x1, c1, c2


def _nearest_step(x, step_size):
    return ((x - 1) // step_size + 1) * step_size


def find_opt_partition(
    stats,
    i: Optional[Union[None, int]] = None,
    j: Optional[Union[None, int]] = None,
    step_size: int = 10,
    delta: Optional[int] = None,
    use_nearest_step: bool = False,
):
    """Find the optimal partition in range (i, j) given the stepsize."""

    # readjust positions to nearest step size to benefit from cache
    if i is None:
        i = 0
    elif use_nearest_step:
        i = _nearest_step(i, step_size)
    if j is None:
        j = len(stats.seq)
    elif use_nearest_step:
        j = min(len(stats.seq), _nearest_step(j, step_size))

    x = np.arange(i, j, step_size, dtype=np.int32)
    _, c1, c2 = get_partitions(stats, x, delta=delta)
    c = c1 + c2
    non_nan_filter = ~np.isnan(c)
    x = x[non_nan_filter]
    c = c[non_nan_filter]

    if not len(c):
        if delta is None:
            return None, None
        else:
            return None, None, None
    cmin = c.min()
    i = np.argmin(c)
    return x[i], cmin


def find_fast_opt_partition(
    stats,
    i=None,
    j=None,
    step_size: int = 100,
    threshold: int = 10,
    delta: Optional[int] = None,
    step_size_refinement: int = 10,
):
    f = functools.partial(find_opt_partition, delta=delta)
    p, pmin = f(stats, i=i, j=j, step_size=step_size)
    if p is None:
        return None, None
    if pmin < threshold:
        return p, pmin
    p2, p2min = f(
        stats,
        max(0, p - step_size),
        min(p + step_size, len(stats.seq)),
        step_size=int(step_size / step_size_refinement),
    )
    if p2 is None:
        return p, pmin
    else:
        return p2, p2min


def find_best_partitions(
    stats: DNAStats,
    threshold: int,
    i=None,
    j=None,
    step_size: int = 100,
    delta: Optional[int] = None,
    partitions=None,
):
    if partitions is None:
        partitions = []
    c = stats.cost(i, j)
    if c < threshold:
        return partitions
    p, pmin = find_fast_opt_partition(stats, i=i, j=j, step_size=step_size, delta=delta)
    if p is None or p in partitions:
        return partitions
    else:
        partitions.append(p)
    c1 = stats.cost(i, p)
    c2 = stats.cost(p, j)
    if c1 > threshold:
        find_best_partitions(
            stats, threshold=threshold, i=i, j=p, partitions=partitions
        )
    if c2 > threshold:
        find_best_partitions(
            stats, threshold=threshold, i=p, j=j, partitions=partitions
        )

    return partitions


def _shift_indices(indices: List[int], origin: int, length: int):
    shifted_indices = []
    for i in indices:
        i += i + origin
        if i >= length:
            i -= length
        shifted_indices.append(i)
    return shifted_indices


def find_by_partitions_for_sequence(
    stats: DNAStats,
    cyclic: bool,
    threshold: int,
    step_size: int = 100,
    delta: Optional[int] = None,
):
    """Approximates the best partitions for a sequence. If cyclic=True, then
    will approximate partitions by also rotating the origin.

    :param stats: DNAStats instance
    :param cyclic: whether the sequence is cyclic
    :param threshold: threshold cost to find a partition
    :param step_size: step size to find partition.
    :return:
    """
    f = functools.partial(
        find_best_partitions, threshold=threshold, step_size=step_size, delta=delta
    )
    partitions = f(stats)
    if cyclic:
        origin = int(len(stats.seq) / 2.0)
        seq = stats.seq
        stats2 = stats.copy_with_new_seq(seq[origin:] + seq[:origin])
        partitions2 = f(stats2)
        partitions += _shift_indices(partitions2, origin, len(seq))
    partitions = sorted(set(partitions))
    return partitions
