from dasi.utils import Span
from dasi import Alignment, Constants, AlignmentException
import pytest


def test_init():
    query_region = Span(0, 1000, 10000)
    subject_region = Span(1000, 2000, 5000)
    alignment = Alignment(query_region, subject_region, Constants.PCR_PRODUCT, "query_key", "subject_key")

    assert alignment.query_key == 'query_key'
    assert alignment.subject_key == 'subject_key'
    assert alignment.query_region == query_region
    assert alignment.subject_region == subject_region


def test_init_raises_alignment_exception():
    """Expect an error if query_region and subject_region have different lengths"""
    query_region = Span(0, 1000, 10000)
    subject_region = Span(1001, 2000, 5000)
    with pytest.raises(AlignmentException):
        Alignment(query_region, subject_region, Constants.PCR_PRODUCT, "query_key", "subject_key")


def test_subregion():
    query_region = Span(0, 1000, 10000)
    subject_region = Span(1000, 2000, 5000)
    alignment = Alignment(query_region, subject_region, Constants.PCR_PRODUCT, "query_key", "subject_key")

    new_alignment = alignment.sub_region(500, 800)

    assert new_alignment.query_key == 'query_key'
    assert new_alignment.subject_key == 'subject_key'
    assert new_alignment.query_region.a == 500
    assert new_alignment.query_region.b == 800
    assert new_alignment.subject_region.a == 1500
    assert new_alignment.subject_region.b == 1800

