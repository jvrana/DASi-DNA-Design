import numpy as np

from dasi.cost.utils import lexargmin


def test_lexargmin():
    a = np.ones((50, 51, 52, 53))  # (20*30*40).reshape(20,30,40)

    b = np.random.randint(0, 1000, size=(50, 51, 52, 53))

    i = lexargmin((a, b), axis=0)
    print(i)
