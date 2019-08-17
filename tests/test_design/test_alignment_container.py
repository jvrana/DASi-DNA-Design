from dasi import Alignment, AlignmentContainer, Constants, AlignmentContainerException
from dasi.utils import Span
import pytest
import random
from uuid import uuid4
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
    return Span(start, end, context_len, cyclic=False)


def random_span_with_len(context_len, l):
    start = random.randint(0, context_len - l)
    end = start + l
    return Span(start, end, context_len)


def random_alignment(type, query_key=None, subject_key=None):
    if query_key is None:
        query_key = str(uuid4())
    if subject_key is None:
        subject_key = str(uuid4())
    query_region = random_span(10000)
    subject_region = random_span_with_len(10000, len(query_region))
    return Alignment(query_region, subject_region, type, query_key, subject_key)


def random_container(num_records, num_alignments, type):
    records = [random_record(10000) for _ in range(num_records)]
    seqdb = {}
    alignments = []
    for r in records:
        seqdb[r.id] = r
    for _ in range(num_alignments):
        j = random.randint(1, len(records)-1)

        align = random_alignment(type, records[0].id, records[j].id)
        alignments.append(align)
    return AlignmentContainer(seqdb, alignments=alignments)


def test_init_raise_error():
    container = random_container(100, 300, Constants.PCR_PRODUCT)
    with pytest.raises(AlignmentContainerException):
        container.alignments = container.alignments + [random_alignment(Constants.PCR_PRODUCT)]


def test_init():
    random_container(100, 300, Constants.PCR_PRODUCT)


def test_different_queries_raises_error():
    raise NotImplementedError


def test_group():
    raise NotImplementedError


def test_expand():
    raise NotImplementedError


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
