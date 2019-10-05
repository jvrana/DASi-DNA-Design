"""Test NumpyDataFrame."""
import numpy as np
import pytest

from dasi.utils import NumpyDataFrame
from dasi.utils import NumpyDataFrameException


@pytest.fixture(scope="function")
def a():
    return NumpyDataFrame({"A": np.arange(10), "B": np.arange(10, 20)})


def test_init():
    a = NumpyDataFrame({"A": np.arange(10), "B": np.arange(10)})
    assert a


def test_empty_init():
    a = NumpyDataFrame()
    assert a.data == {}


def test_init_raises():
    with pytest.raises(NumpyDataFrameException):
        NumpyDataFrame({"A": np.arange(10), "B": np.arange(9)})


def test_set_item(a):
    a[2] = 100.0
    assert a.data["B"][2] == 100.0
    assert np.all(a.data["A"] == np.array([0, 1, 100.0, 3, 4, 5, 6, 7, 8, 9]))


@pytest.mark.parametrize("col", ["A", "B", "C"])
def test_new_col(a, col):
    a.col[col] = np.arange(100, 110)
    assert np.all(a.data[col] == np.arange(100, 110))


def test_print(a):
    print(a)


def test_shape(a):
    assert a.shape == (10,)


def test_reshape(a):
    b = a.reshape((10, 1))
    assert a is not b
    assert a.shape == (10,)
    assert b.shape == (10, 1)


def test_squeeze(a):
    a = a.reshape((10, 1))
    assert a.shape == (10, 1)
    assert a.apply(np.squeeze).shape == (10,)


def test_apply_with_kwargs(a):
    df = NumpyDataFrame({"A": np.arange(10), "B": np.arange(10, 20)})
    df = df.reshape((-1, 1))

    df1 = df.apply(np.sum, axis=1)
    assert np.all(df1.data["A"] == np.arange(10))
    assert np.all(df1.data["B"] == np.arange(10, 20))

    df2 = df.apply(np.sum, axis=0)
    assert np.all(df2.data["A"] == np.arange(10).sum(axis=0))
    assert np.all(df2.data["B"] == np.arange(10, 20).sum(axis=0))


def test_len(a):
    assert len(a) == 10
    b = a.group_apply([a, a], np.hstack)
    print(b.shape)
    assert len(b) == 20


def test_copy(a):
    b = a.copy()
    for k in a.columns:
        assert a.data[k] is not b.data[k]


def test_concat(a):
    c = a.copy()
    d = NumpyDataFrame.concat([a, c])
    assert d.shape == (20,)


def test_concat_raises(a):
    c = a.copy()
    c.col["C"] = np.arange(100, 110)
    with pytest.raises(NumpyDataFrameException):
        d = NumpyDataFrame.concat([a, c])


def test_concat_fills_missing(a):
    c = a.copy()
    c.col["C"] = np.arange(100, 110)
    d = NumpyDataFrame.concat([a, c], fill_value=np.inf)
    assert np.all(d.data["C"] == np.array([np.inf] * 10 + list(range(100, 110))))


class TestUpdateMerge:
    def test_update(self, a):
        b = NumpyDataFrame({"C": np.arange(30, 40)})
        a.update(b)
        assert "C" in a.columns
        print(a)
        assert np.all(a.data["A"] == np.arange(10))
        assert np.all(a.data["C"] == np.arange(30, 40))
        assert np.all(a.data["B"] == np.arange(10, 20))

    def test_update_raises(self, a):
        b = NumpyDataFrame({"C": np.arange(30, 41)})
        with pytest.raises(NumpyDataFrameException):
            a.update(b)

    def test_merge(self, a):
        b = NumpyDataFrame({"C": np.arange(30, 40)})
        c = NumpyDataFrame.merge((a, b))
        assert "C" in c.columns
        assert np.all(c.data["A"] == np.arange(10))
        assert np.all(c.data["C"] == np.arange(30, 40))
        assert np.all(c.data["B"] == np.arange(10, 20))


class Slicer:
    def __getitem__(self, item):
        return item


slicer = Slicer()


class TestIndexing:
    @pytest.mark.parametrize(
        "s",
        [
            slicer[0],
            slicer[0:10],
            slicer[5:],
            slicer[:6],
            slicer[:-1],
            slicer[1:],
            slicer[-2],
        ],
    )
    def test_index(self, a, s):
        df = a[s]
        assert np.all(df.data["A"] == np.arange(10)[s])
        assert np.all(df.data["B"] == np.arange(10, 20)[s])

    @pytest.mark.parametrize(
        "s",
        [
            slicer[0],
            slicer[0:10],
            slicer[5:],
            slicer[:6],
            slicer[:-1],
            slicer[1:],
            slicer[-2],
        ],
    )
    def test_2D_index(self, a, s):
        a = a.reshape((10, 1))
        df = a[s]
        assert np.all(df.data["A"] == np.arange(10).reshape(10, 1)[s])
        assert np.all(df.data["B"] == np.arange(10, 20).reshape(10, 1)[s])

    def test_col_index(self, a):
        assert a.col["A"].columns == ("A",)
        assert a.col["B"].columns == ("B",)
        assert a.col["A", "B"].columns == ("A", "B")
        assert a.col["B", "A"].columns == ("A", "B")
        assert a.col["A"]

        assert np.all(a.col["A", "B"].data["A"] == np.arange(10))
        assert np.all(a.col["A", "B"].data["B"] == np.arange(10, 20))
        assert "A" in a.col
        assert "B" in a.col
        assert "C" not in a.col
        assert "B" not in a.col["A"].col


class TestOperations:
    def test_add_int(self, a):
        df = a + 10
        assert np.all(df.data["A"] == np.arange(10) + 10)

    def test_add_df(self, a):
        df = a + a
        assert np.all(df.data["A"] == np.arange(10) + np.arange(10))

    def test_mul_float(self, a):
        df = a * 1.5
        assert np.all(df.data["A"] == np.arange(10) * 1.5)

    def test_mul_df(self, a):
        df = a * a
        assert np.all(df.data["A"] == np.arange(10) * np.arange(10))

    def test_pow_float(self, a):
        df = a ** 1.5
        assert np.all(df.data["A"] == np.arange(10) ** 1.5)

    def test_pow_df(self, a):
        df = a ** a
        assert np.all(df.data["A"] == np.arange(10) ** np.arange(10))

    def test_div_float(self, a):
        df = a ** 1.5
        assert np.all(df.data["A"] == np.arange(10) ** 1.5)

    def test_div_df(self, a):
        df = a / (a + 1)
        assert np.all(df.data["A"] == np.arange(10) / (np.arange(10) + 1))

    def test_neg(self, a):
        df = -a
        assert np.all(df.data["A"] == -np.arange(10))

    def test_sub_int(self, a):
        df = a - 10
        assert np.all(df.data["A"] == np.arange(10) - 10)

    def test_sub_df(self, a):
        df = a - a
        assert np.all(df.data["A"] == np.arange(10) - np.arange(10))


@pytest.mark.parametrize("shape", [tuple(), (1,), (5,), (1, 5), (5, 1), (5, 5)])
def test_can_slice(shape):
    a = NumpyDataFrame({"A": np.ones(shape), "B": np.zeros(shape)})
    if len(shape) == 0:
        with pytest.raises(IndexError):
            a[0]
    else:
        assert a[0] is not None


@pytest.mark.parametrize("shape", [tuple(), (1,), (5,), (1, 5), (5, 1), (5, 5)])
def test_str(shape):
    a = NumpyDataFrame({"A": np.ones(shape), "B": np.zeros(shape)})
    print(str(a))


@pytest.mark.parametrize("shape", [tuple(), (1,), (5,), (1, 5), (5, 1), (5, 5)])
def test_repr(shape):
    a = NumpyDataFrame({"A": np.ones(shape), "B": np.ones(shape)})
    print(a.__repr__())


@pytest.mark.parametrize("shape", [tuple(), (1,), (5,), (1, 5), (5, 1)])
def test_to_df(shape):
    a = NumpyDataFrame({"A": np.ones(shape), "B": np.zeros(shape)})
    assert a.shape == shape
    print(a.to_df())


@pytest.mark.parametrize("shape", [(2, 5), (5, 2)])
def test_to_df_raises(shape):
    a = NumpyDataFrame({"A": np.ones(shape), "B": np.zeros(shape)})
    with pytest.raises(NumpyDataFrameException):
        a.to_df()
