import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dasi.utils.region import Span


class TestInit:
    @pytest.mark.parametrize("a", list(range(10)))
    @pytest.mark.parametrize("b", list(range(10)))
    @pytest.mark.parametrize("cyclic", [True, False])
    def test_all_init(self, a, b, cyclic):
        length = 20
        if a > b and not cyclic:
            with pytest.raises(IndexError):
                Span(a, b, length, cyclic)
        else:
            s = Span(a, b, length, cyclic)
            assert s.a == a
            assert s.b == b
            assert s.context_length == length
            assert s.context_length == length
            if cyclic and a > b:
                assert len(s) == length - a + b
            else:
                assert len(s) == b - a

    def test_init_should_raise(self):
        with pytest.raises(IndexError):
            Span(10, 10, 10, False)

    def test_init_should_raise2(self):
        with pytest.raises(IndexError):
            Span(9408, 4219, 9408, True)

    def test_init_linear(self):
        s = Span(10, 80, 100, False)
        assert s.a == 10
        assert s.b == 80
        assert s.context_length == 100
        assert s.cyclic is False

    def test_init_cyclic(self):
        s = Span(10, 5, 100, True)
        assert s.a == 10
        assert s.b == 5
        assert s.context_length == 100
        assert s.cyclic

    def test_init_linear_raises(self):
        Span(10, 5, 100, True)
        with pytest.raises(IndexError):
            Span(10, 5, 100, False)

    def test_invalid_cyclic(self):
        print(Span(0, 10, 10, cyclic=True))
        with pytest.raises(IndexError):
            print(Span(0, 10, 9, cyclic=True, strict=True))
        print(Span(0, 10, 9, True, strict=False))


def test_len_linear():
    s = Span(5, 20, 100, False)
    assert len(s) == 15


def test_len_cyclic():
    s = Span(90, 10, 100, True)
    assert len(s) == 20


def test_iter_linear():
    s = Span(10, 20, 100, False)
    assert list(s) == list(range(10, 20))


def test_iter_cyclic():
    s = Span(90, 10, 100, True)
    assert list(s) == list(range(90, 100)) + list(range(10))


def test_str():
    print(Span(90, 10, 100, True))


def test_eq():
    s1 = Span(10, 90, 100, True)
    s2 = Span(10, 90, 100, True)
    assert s1 == s2

    s1 = Span(10, 91, 100, True)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(11, 90, 100, True)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(10, 90, 100, False)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(10, 90, 100, True)
    s2 = Span(10, 90, 101, True)
    assert not s1 == s2


class TestContains:
    @pytest.mark.parametrize("i", range(-20, 20))
    def test_contains_index(self, i):
        s = Span(10, 50, 100, True)
        if 10 <= i < 50:
            assert i in s
        else:
            assert not i in s

    class TestContainsCyclic:
        @pytest.mark.parametrize("i", range(90, 100))
        def test_contains_cyclic_1(self, i):
            s = Span(90, 10, 100, True)
            assert i in s

        @pytest.mark.parametrize("i", range(0, 10))
        def test_contains_cyclic_2(self, i):
            s = Span(90, 10, 100, True)
            assert i in s

        @pytest.mark.parametrize("i", range(85, 90))
        def test_not_contains_cyclic_1(self, i):
            s = Span(90, 10, 100, True)
            assert i not in s

        @pytest.mark.parametrize("i", range(11, 20))
        def test_not_contains_cyclic_2(self, i):
            s = Span(90, 10, 100, True)
            assert i not in s

    class TestContainsSpan:
        def test_contains_region(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(20, 30, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_empty(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(30, 30, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_empty2(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(10, 10, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_self(self):
            s1 = Span(10, 50, 100, True)
            assert s1 in s1

        def test_does_not_contain(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(11, 50, 100, True)
            assert s1 not in s2
            assert s2 in s1

            s1 = Span(10, 50, 100, True)
            s2 = Span(10, 49, 100, True)
            assert s1 not in s2
            assert s2 in s1

        def test_cyclic_contains(self):
            s1 = Span(80, 20, 100, True)
            s2 = Span(85, 10, 100, True)
            assert s2 in s1
            assert s1 not in s2

        def test_contains_example(self):
            r1 = Span(5947, 4219, 9408, True)
            with pytest.raises(IndexError):
                r2 = Span(9408, 4219, 9408, True)

        def test_contains_example2(self):
            r1 = Span(59, 42, 94, True)
            with pytest.raises(IndexError):
                r2 = Span(94, 42, 94, True)


class TestIntersection:
    def x(self, a1, b1, a2, b2):
        s1 = Span(a1, b1, 100, True)
        s2 = Span(a2, b2, 100, True)
        return s1, s2

    def test_intersection(self):
        s1 = Span(10, 50, 100, True)
        s2 = Span(40, 60, 100, True)

        sliced = s1.intersection(s2)
        assert sliced.a == 40
        assert sliced.b == 50

        sliced = s2.intersection(s1)
        assert sliced.a == 40
        assert sliced.b == 50

    @pytest.mark.parametrize("i", [-1, 0, 1, 2])
    def test_no_intersection(self, i):
        s1, s2 = self.x(10, 50, 50 + i, 80)
        if i >= 0:
            assert not s1.intersection(s2)
            assert not s2.intersection(s1)
        else:
            assert s1.intersection(s2)
            assert s2.intersection(s1)

    def test_encompassed_intersection(self):
        s1, s2 = self.x(10, 50, 20, 30)
        sliced = s1.intersection(s2)
        assert sliced.a == 20
        assert sliced.b == 30

        sliced = s2.intersection(s1)
        assert sliced.a == 20
        assert sliced.b == 30


class TestSlice:
    @pytest.mark.parametrize("i", list(range(-20, 20)))
    @pytest.mark.parametrize("cyclic", [True, False])
    def test_indexing(self, i, cyclic):

        s = Span(5, 15, 20, cyclic)
        assert len(s) == 10
        if (i < -len(s) or i >= len(s)) and not cyclic:
            with pytest.raises(IndexError):
                print(s[i])
        else:
            if i < 0:
                assert s[i] == (i + 15) % 20
            else:
                assert s[i] == (i + 5) % 20

    @pytest.mark.parametrize("i", list(range(10, 50, 5)))
    @pytest.mark.parametrize("j", list(range(60, 80, 5)))
    @pytest.mark.parametrize("cyclic", [True, False])
    def test_valid_slices(self, i, j, cyclic):
        s = Span(10, 100, 200, cyclic)
        sliced = s[i:j]
        assert len(sliced) == j - i

    @pytest.mark.parametrize("i", range(10))
    def test_slice_and_invert_slice_should_total_length(self, i):
        s = Span(10, 100, 200, True)
        if i != 4:
            x1 = s[4:i]
            x2 = s[i:4]
            assert len(x1) + len(x2) == 200

    def test_open_ended_slices(self):
        s = Span(0, 10, 20, cyclic=False)
        assert list(s[:10]) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert list(s[:-10]) == []
        with pytest.raises(IndexError):
            assert list(s[10:]) == []
        assert list(s[-10:]) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    def test_open_ended_slices_cyclic(self):
        s = Span(0, 10, 20, cyclic=True)
        assert list(s[:10]) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert list(s[:-10]) == []
        assert list(s[10:]) == []
        assert list(s[-10:]) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    @pytest.mark.parametrize("i", list(range(-20, 20)))
    def test_open_ended_slice_left(self, i):
        s = Span(10, 20, 200, False)
        if i > 10 or i < -10:
            with pytest.raises(IndexError):
                s[:i]
        else:
            sliced = s[:i]
            assert sliced.a == 10
            if i < 0:
                assert sliced.b == 20 + i
            else:
                assert sliced.b == 10 + i

    @pytest.mark.parametrize("i", list(range(-20, 20)))
    def test_open_ended_slice_right(self, i):
        s = Span(10, 20, 200, False)
        if i >= 10 or i < -10:
            with pytest.raises(IndexError):
                assert not s[i:]
        else:
            sliced = s[i:]
            assert sliced.b == 20
            if i < 0:
                assert sliced.a == 20 + i
            else:
                assert sliced.a == 10 + i

    def test_invalid_slice(self):
        s = Span(90, 10, 100, True)
        sliced = s[:15]
        assert sliced.a == 90
        assert sliced.b == 5

        s = Span(90, 10, 100, True)
        sliced = s[-15:]
        assert sliced.a == 95
        assert sliced.b == 10

    def test_copying(self):
        s = Span(5, 15, 20, True)
        sliced = s[:]
        assert sliced.a == s.a
        assert sliced.b == s.b

    def test_copy_new_ab(self):
        s = Span(5, 15, 20, True)
        copied = s[2, 16]
        assert copied.a == 2
        assert copied.b == 16

    def test_invert_cyclic(self):

        s = Span(5, 15, 20, True)
        i1 = s.invert()[0]

        assert i1.a == 15
        assert i1.b == 5

    def test_invert_linear(self):
        s = Span(5, 15, 20, False)
        s1, s2 = s.invert()
        assert s1.a == 0
        assert s1.b == 5
        assert s2.a == 15
        assert s2.b == 20


class TestDifference:
    def test_linear_diff(self):

        s1 = Span(20, 80, 100, True)
        s2 = Span(30, 50, 100, True)

        diff = s1.differences(s2)
        assert len(diff) == 2
        assert diff[0].a == 20
        assert diff[0].b == 30
        assert diff[1].a == 50
        assert diff[1].b == 80

    def test_overhang_left_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(10, 30, 100, True)

        diff = s1.differences(s2)
        assert len(diff) == 1
        assert diff[0].a == 30
        assert diff[0].b == 80

        diff = s2.differences(s1)
        print(diff)
        assert len(diff) == 1
        assert diff[0].a == 10
        assert diff[0].b == 20

    def test_overhang_right_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(70, 90, 100, True)
        diff = s1.differences(s2)
        assert len(diff) == 1
        assert diff[0].a == 20
        assert diff[0].b == 70

        diff = s2.differences(s1)
        print(diff)
        assert len(diff) == 1
        assert diff[0].a == 80
        assert diff[0].b == 90

    def test_no_overlap_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(90, 95, 100, True)
        diff = s1.differences(s2)
        assert len(diff) == 1
        assert s1.a == 20
        assert s1.b == 80


class TestConsecutive:
    def test_consecutive(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(81, 95, 100, True)
        assert not s1.consecutive(s2)
        assert not s2.consecutive(s1)

    def test_not_consecutive(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(80, 95, 100, True)
        assert s1.consecutive(s2)
        assert not s2.consecutive(s1)

    def test_consecutive_over_origin1(self):
        s1 = Span(7, 0, 10, True)
        s2 = Span(0, 2, 10, True)
        assert s1.consecutive(s2)
        assert not s2.consecutive(s1)


def test_connecting_span():
    s1 = Span(20, 50, 100, False)
    s2 = Span(60, 75, 100, False)
    s3 = s1.connecting_span(s2)
    assert s3.a == 50
    assert s3.b == 60


def test_connecting_span_cyclic():
    s1 = Span(10, 20, 100, True)
    s2 = Span(80, 90, 100, True)
    s3 = s1.connecting_span(s2)
    assert s3.a == 20
    assert s3.b == 80

    s4 = s2.connecting_span(s1)
    assert s4.a == 90
    assert s4.b == 10


def test_connecting_span_at_origin():

    s1 = Span(50, 60, 100, True)
    s2 = Span(0, 30, 100, True)
    s3 = s1.connecting_span(s2)
    assert s3.a == 60
    assert s3.b == 0


def test_connectin_span_over_origin():

    s1 = Span(50, 60, 100, True)
    s2 = Span(5, 30, 100, True)
    s3 = s1.connecting_span(s2)
    assert s3.a == 60
    assert s3.b == 5


def test_self_connecting_span():
    s1 = Span(50, 60, 100, True)
    s2 = s1.connecting_span(s1)
    assert s2.a == 60
    assert s2.b == 50


def test_connectin_span_with_overlap():

    s1 = Span(10, 30, 100, True)
    s2 = Span(20, 50, 100, True)
    assert not s1.connecting_span(s2)


def test_connectin_span_consecurive_is_empty():

    s1 = Span(10, 30, 100, True)
    s2 = Span(30, 50, 100, True)
    assert not s1.overlaps_with(s2)
    s3 = s1.connecting_span(s2)
    assert len(s3) == 0


def test_connecting_span_linear_no_span():
    s1 = Span(10, 99, 100, False)
    s2 = Span(0, 10, 100, False)
    s3 = s1.connecting_span(s2)
    assert s3 is None


def test_invalid_span():
    s1 = Span(5947, 4219, 10000, cyclic=True)
    with pytest.raises(IndexError):
        s1.sub(28, 5989)


def test_ranges_ignore_wraps():

    s1 = Span(1, 5, 10, cyclic=True)
    assert s1.ranges() == [(1, 5)]
    assert s1.ranges(ignore_wraps=True) == [(1, 5)]

    s1 = Span(1, 15, 10, cyclic=True)
    assert s1.ranges() == [(1, 10), (0, 5)]
    assert s1.ranges(ignore_wraps=True) == [(1, 5)]

    s1 = Span(1, 25, 10, cyclic=True)
    assert s1.ranges() == [(1, 10), (0, 10), (0, 5)]
    assert s1.ranges(ignore_wraps=True) == [(1, 5)]

    s2 = Span(8, 2, 10, cyclic=True)
    assert s2.ranges() == [(8, 10), (0, 2)]
    assert s2.ranges(ignore_wraps=True) == [(8, 10), (0, 2)]

    # TODO: this is a strange behavior?
    s2 = Span(8, 12, 10, cyclic=True)
    assert s2.ranges() == [(8, 10), (0, 10), (0, 2)]
    assert s2.ranges(ignore_wraps=True) == [(8, 10), (0, 2)]

    s2 = Span(8, 22, 10, cyclic=True)
    assert s2.ranges() == [(8, 10), (0, 10), (0, 10), (0, 2)]
    assert s2.ranges(ignore_wraps=True) == [(8, 10), (0, 2)]


class TestGetSlice:
    def test_get_slice_str():
        span = Span(8, 22, 10, cyclic=True)
        x = "abcdefghijklmnopqrstuv"
        y = span.get_slice(x)

    def test_get_slice_list_of_ints(self):
        span = Span(8, 22, 10, cyclic=True)
        x = list(range(20, 140))
        y = span.get_slice(x)
        print(list(span.get_slice_iter(x)))
        print(y)
        assert y == [x[_s] for _s in list(span)]

        span = Span(8, 22, 10, cyclic=True)
        x = tuple(range(20, 140))
        y = span.get_slice(x)
        assert y == tuple([x[_s] for _s in list(span)])

    def test_get_slice_seqrecord(self):
        seq = SeqRecord(Seq("AGGGTGTGTGCGA"))
        span = Span(8, 22, len(seq), cyclic=True)

        seq2 = span.get_slice(seq)
        print(seq2)
