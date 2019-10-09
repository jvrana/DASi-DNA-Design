import pytest

from dasi.alignments import Alignment
from dasi.constants import Constants
from dasi.exceptions import AlignmentException
from dasi.utils.region import Region


def test_init():
    query_region = Region(0, 1000, 10000)
    subject_region = Region(1000, 2000, 5000)
    alignment = Alignment(
        query_region, subject_region, Constants.PCR_PRODUCT, "query_key", "subject_key"
    )

    assert alignment.query_key == "query_key"
    assert alignment.subject_key == "subject_key"
    assert alignment.query_region == query_region
    assert alignment.subject_region == subject_region


def test_init_raises_alignment_exception():
    """Expect an error if query_region and subject_region have different
    lengths."""
    query_region = Region(0, 1000, 10000)
    subject_region = Region(1001, 2000, 5000)
    with pytest.raises(AlignmentException):
        Alignment(
            query_region,
            subject_region,
            Constants.PCR_PRODUCT,
            "query_key",
            "subject_key",
        )


class TestSubRegion:
    def test_subregion(self):
        query_region = Region(0, 1000, 10000)
        subject_region = Region(1000, 2000, 5000)
        alignment = Alignment(
            query_region,
            subject_region,
            Constants.PCR_PRODUCT,
            "query_key",
            "subject_key",
        )

        new_alignment = alignment.sub_region(500, 800)

        assert new_alignment.query_key == "query_key"
        assert new_alignment.subject_key == "subject_key"
        assert new_alignment.query_region.a == 500
        assert new_alignment.query_region.b == 800
        assert new_alignment.subject_region.a == 1500
        assert new_alignment.subject_region.b == 1800

    @pytest.mark.parametrize(
        "x",
        [
            (100, 1000, 10000, 1100, 2000, 5000, 100, -100),
            (100, 1000, 10000, 1100, 2000, 5000, 200, -200),
            (9000, 1000, 10000, 1000, 3000, 5000, 200, -200),
            (100, 2000, 10000, 4100, 1000, 5000, 100, -100),
            (9000, 1000, 10000, 4000, 1000, 5000, 200, -200),
        ],
        ids=[
            "basic1",
            "basic2",
            "query_over_origin",
            "subject_over_origin",
            "both_over_origin",
        ],
    )
    def test_subregion2(self, x):
        query_region = Region(x[0], x[1], x[2], cyclic=True)
        subject_region = Region(x[3], x[4], x[5], cyclic=True)
        alignment = Alignment(
            query_region,
            subject_region,
            Constants.PCR_PRODUCT,
            "query_key",
            "subject_key",
        )

        _q = query_region[x[6] : x[7]]
        _s = subject_region[x[6] : x[7]]
        new_alignment = alignment.sub_region(_q.a, _q.b)

        assert new_alignment.query_key == "query_key"
        assert new_alignment.subject_key == "subject_key"
        assert new_alignment.query_region.a == _q.a
        assert new_alignment.query_region.b == _q.b
        assert new_alignment.subject_region.a == _s.a
        assert new_alignment.subject_region.b == _s.b

    @pytest.mark.parametrize(
        "x",
        [
            (100, 1000, 10000, 1100, 2000, 5000, 100, -100),
            (100, 1000, 10000, 1100, 2000, 5000, 200, -200),
            (9000, 1000, 10000, 1000, 3000, 5000, 200, -200),
            (100, 2000, 10000, 4100, 1000, 5000, 100, -100),
            (9000, 1000, 10000, 4000, 1000, 5000, 200, -200),
        ],
        ids=[
            "basic1",
            "basic2",
            "query_over_origin",
            "subject_over_origin",
            "both_over_origin",
        ],
    )
    def test_subregion_reversed(self, x):
        query_region = Region(x[0], x[1], x[2], cyclic=True)
        subject_region = Region(x[3], x[4], x[5], cyclic=True, direction=-1)
        alignment = Alignment(
            query_region,
            subject_region,
            Constants.PCR_PRODUCT,
            "query_key",
            "subject_key",
        )

        _q = query_region[x[6] : x[7]]
        _s = subject_region[-x[7] : -x[6]]
        sub = alignment.sub_region(_q.a, _q.b)

        assert sub.subject_region.a == _s.a
        assert sub.subject_region.b == _s.b
