from dasi.utils import SpanDictionary
import pytest


def test_init_span_dictionary():
    s = SpanDictionary()
    assert not s
    s[1, 2] = 1
    assert s


def test_to_dict1():
    s = SpanDictionary({(1, 2): 1})

    assert s.to_dict() == {(1, 2): 1}


def test_to_dict2():
    s = SpanDictionary({(1, 2): 1, (2, 5): 5})

    assert s.to_dict() == {(1, 2): 1, (2, 5): 5}


def test_get_item():
    s = SpanDictionary({(4, 7): 10})

    for i in range(10):
        if 4 <= i < 7:
            assert s[i] == [10]
        else:
            assert not s[i]


@pytest.mark.parametrize("i", list(range(10)))
def test_overlapping_spans(i):
    s = SpanDictionary({(4, 7): 10, (5, 6): 15})

    if 4 <= i < 7:
        if i == 5:
            assert s[i] == [10, 15]
        else:
            assert s[i] == [10]
    else:
        assert not s[i]


def test_iter():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})

    assert list(s) == [4, 5, 6, 9, 10]


def test_items():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})

    assert list(s.items()) == [
        (4, [10]),
        (5, [10, 15]),
        (6, [10]),
        (9, [40]),
        (10, [40]),
    ]


def test_min():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})

    assert s.min() == 4


def test_max():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})

    assert s.max() == 11


def test_ranges():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})
    assert list(s.ranges()) == [(4, 7), (5, 6), (9, 11)]


def test_len():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})
    assert len(s) == 3


def test_keys():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})
    assert list(s.keys()) == [(4, 7), (5, 6), (9, 11)]


def test_values():
    s = SpanDictionary({(4, 7): 10, (5, 6): 15, (9, 11): 40})
    assert list(s.values()) == [10, 15, 40]
