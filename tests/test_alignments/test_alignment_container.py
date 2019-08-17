from dasi import Alignment, AlignmentContainer, Constants, AlignmentContainerException
from dasi.utils import Region
import pytest
import random
from uuid import uuid4
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from copy import copy



def random_seq(len):
    bases = "AGTC"

    i = random.randint(0, 3)
    seq = ''
    for _ in range(len):
        seq += bases[i]
    return seq


def random_record(len):
    return SeqRecord(Seq(random_seq(len)), id=str(uuid4()))\


def random_span(context_len):
    random_len = random.randint(100, context_len - 100)
    start = random.randint(0, context_len - random_len)
    end = start + random_len
    return Region(start, end, context_len, cyclic=False, direction=1)


def random_span_with_len(context_len, l):
    start = random.randint(0, context_len - l)
    end = start + l
    return Region(start, end, context_len, direction=1)


def random_alignment(type, query_key=None, subject_key=None, span_len=None, context_len=10000):
    if query_key is None:
        query_key = str(uuid4())
    if subject_key is None:
        subject_key = str(uuid4())
    query_region = random_span(context_len)
    if span_len is None:
        span_len = len(query_region)
    subject_region = random_span_with_len(context_len, span_len)
    return Alignment(query_region, subject_region, type, query_key, subject_key)


def random_container(num_records, num_alignments, type, context_len=10000, span_len=None):
    records = [random_record(context_len) for _ in range(num_records)]
    seqdb = {}
    alignments = []
    for r in records:
        seqdb[r.id] = r
    for _ in range(num_alignments):
        j = random.randint(1, len(records)-1)

        align = random_alignment(type, records[0].id, records[j].id, context_len=context_len, span_len=span_len)
        alignments.append(align)
    return AlignmentContainer(seqdb, alignments=alignments)



def test_init():
    container = random_container(100, 300, Constants.PCR_PRODUCT)
    assert len(container) == 300


def test_init_raise_error_with_different_queries():
    container = random_container(100, 300, Constants.PCR_PRODUCT)
    with pytest.raises(AlignmentContainerException):
        container.alignments = container.alignments + [random_alignment(Constants.PCR_PRODUCT)]


def test_group():
    a1 = random_alignment(Constants.PCR_PRODUCT)
    a2 = copy(a1)
    a3 = random_alignment(Constants.FRAGMENT)

    alignments = [a1, a2, a3]

    groups = AlignmentContainer.group(alignments)
    assert len(groups) == 2


def test_group_by_types():
    """Tests if we correctly gather alignments by their type"""
    container = random_container(10, 3, Constants.PCR_PRODUCT)

    # add two FRAGMENT types
    align = container.alignments[0]
    align.type = Constants.FRAGMENT
    copied = copy(align)
    container.alignments.append(copied)

    # check if we retrieve the correct group
    assert len(container.get_groups_by_types(Constants.FRAGMENT)) == 1
    assert len(container.get_groups_by_types(Constants.PCR_PRODUCT)) >= 1
    assert len(container.get_groups_by_types(Constants.FRAGMENT)[0].alignments) == 2


def test_expand():
    container = random_container(100, 100, Constants.PCR_PRODUCT)
    container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))




def test_expand_pcr_products():
    raise NotImplementedError


def test_expand_primer_pairs():
    raise NotImplementedError


def test_groups_by_types():
    raise NotImplementedError


def test_alignments_by_types():
    raise NotImplementedError


def test_types():
    raise NotImplementedError


def test_annotate_fragments():
    raise NotImplementedError
