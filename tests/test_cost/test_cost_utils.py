from dasi.cost.utils import lexargmin
import numpy as np


def test_lexargmin():

    a = np.ones((50,100,100, 100)) # (20*30*40).reshape(20,30,40)

    b = np.random.randint(0, 1000, size=(50, 100, 100, 100))

    i = lexargmin((a, b), axis=0)
    print(i)