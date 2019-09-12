import numpy as np


def test_benchmark_choose(benchmark):

    A = np.ones((49, 50))
    B = np.random.randint(0, 10, (49, 50))

    benchmark(np.choose, A > B, (A, B))


def test_benchmark_diy(benchmark):
    A = np.ones((49, 50))
    B = np.random.randint(0, 10, (49, 50))

    def func(A, B):
        a = A.ravel()
        b = B.ravel()

        c = a < b
        d = np.ma.masked_array(a, c)

    benchmark(func, A, B)
