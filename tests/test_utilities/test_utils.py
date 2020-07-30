import numpy as np

from dasi.utils import bisect_slice_between
from dasi.utils import chunkify
from dasi.utils import sort_cycle


def test_sort_cycle():

    a = [5, 6, 7, 1, 2, 3]
    b = [1, 2, 3, 5, 6, 7]
    assert sort_cycle(a) == b


def test_sort_cycle2():
    a = "justin"
    b = "injust"
    s = "{} is {}".format(a, sort_cycle(b))
    assert s == "justin is injust"


def test_benchmark_bisect_slice_between(benchmark):
    a = np.random.randint(0, 1000, 1000)
    a.sort()

    def f():
        b = np.random.randint(0, 1000, 2)
        return bisect_slice_between(a, a, b[0], b[1])

    benchmark(f)


def test_benchmark_bisect_slice_between(benchmark):
    a = np.random.randint(0, 1000, 1000)
    a.sort()

    def f():
        b = np.random.randint(0, 1000, 2)
        return bisect_slice_between(a, a, b[0], b[1])

    benchmark(f)


def test_chunkify1():
    x = list(chunkify([1, 2, 3, 4], 1))
    assert x == [[1], [2], [3], [4]]


def test_chunkify2():
    x = list(chunkify([1, 2, 3, 4, 5], 2))
    assert x == [[1, 2], [3, 4], [5]]


def test_chunkify3():
    x = list(chunkify([1, 2, 3, 4, 5], 6))
    assert x == [[1, 2, 3, 4, 5]]
