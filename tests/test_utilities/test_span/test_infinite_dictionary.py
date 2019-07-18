from shoestring.utils import InfiniteDict
import pytest

def test_ceil_false():
    d = InfiniteDict()
    d[1] = 10
    d[100] = 200

    assert d[-1] == 10
    assert d[0] == 10
    assert d[1] == 10
    assert d[2] == 10
    assert d[3] == 10

    assert d[10] == 10
    assert d[99] == 10
    assert d[100] == 200
    assert d[101] == 200
    assert d[102] == 200


@pytest.mark.parametrize('ceil', [True, False])
@pytest.mark.parametrize('i', range(-300, 300, 10))
def test_infinite(ceil, i):
    d = InfiniteDict(ceil=ceil)
    d[1] = 10
    d[100] = 200
    d[i]


def test_ceil_true():
    d = InfiniteDict(ceil=True)
    d[1] = 10
    d[100] = 200

    assert d[-1] == 10
    assert d[0] == 10
    assert d[1] == 10
    assert d[2] == 200
    assert d[3] == 200
    assert d[10] == 200
    assert d[99] == 200
    assert d[100] == 200
    assert d[101] == 200
    assert d[102] == 200

