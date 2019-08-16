from dasi import Alignment, AlignmentContainer, Constants
from dasi.utils import Span
import pytest
import random
from uuid import uuid4


def random_span(context_len):
    random_len = random.randint(100, context_len - 100)
    start = random.randint(0, context_len-random_len)
    end = start + random_len
    return Span(start, end, context_len, cyclic=False)


def random_span_with_len(context_len, l):
    start = random.randint(0, context_len - l)
    end = start + l
    return Span(start, end, context_len)


def random_alignment(type, query_key=None, subject_key=None):
    if query_key is None:
        query_key = str(uuid4())
    if subject_key is None:
        subjet_key = str(uuid4())
    query_region = random_span(10000)
    subject_region = random_span_with_len(10000, len(query_region))
    return Alignment(query_region, subject_region, type, query_key, subject_key)


def test():
    random_alignment(Constants.PCR_PRODUCT)
