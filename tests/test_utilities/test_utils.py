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
