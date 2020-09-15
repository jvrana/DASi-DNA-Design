import pytest

from dasi.utils.region import Span


@pytest.mark.parametrize("index", range(3))
def test_index_span(index):
    s = Span(90, 10, 100, True, index)
    assert len(s) == 20


@pytest.mark.parametrize("index", range(4))
def test_t(index):
    s = Span(99, 10, 100, True, index)
    assert s.t(-1) == 100 + index - 1
    assert s.t(0) == index
    assert s.t(1) == index + 1


def test_list():
    s = Span(99, 10, 100, True, 0)
    assert list(s) == [99, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    s = Span(99, 10, 100, True, 1)
    assert list(s) == [99, 100, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    s = Span(99, 10, 100, True, 2)
    assert list(s) == [99, 100, 101, 2, 3, 4, 5, 6, 7, 8, 9]


@pytest.mark.parametrize("index", range(5))
def test_indexing(index):
    s = Span(99, 10, 100, True, index)
    assert s[0] == 99
    assert s[-1] == 9
    assert s[-2] == 8


@pytest.mark.parametrize("length", range(1, 5))
@pytest.mark.parametrize("index", range(-5, 5))
def test_slicing_open_left2(length, index):
    s = Span(90, 10, 100, True, index)
    print(s[-length:])
    assert len(s[-length:]) == length


def test_full():
    s = Span(1, 9410, 9409, False, index=1)
    assert len(s) == 9409
