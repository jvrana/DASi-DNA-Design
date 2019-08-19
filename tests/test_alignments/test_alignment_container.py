from dasi.alignments import Alignment, AlignmentContainer
from dasi.constants import Constants
from dasi.exceptions import AlignmentContainerException
from dasi.utils import Region
import pytest
import random
from uuid import uuid4
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from copy import copy, deepcopy


def random_seq(len):
    bases = "AGTC"

    i = random.randint(0, 3)
    seq = ''
    for _ in range(len):
        seq += bases[i]
    return seq


def random_record(len):
    return SeqRecord(Seq(random_seq(len)), id=str(uuid4()))


def random_span(context_len):
    random_len = random.randint(100, context_len - 100)
    start = random.randint(0, context_len - random_len)
    end = start + random_len
    return Region(start, end, context_len, cyclic=True, direction=1)


def random_span_with_len(context_len, l):
    start = random.randint(0, context_len - l)
    end = start + l
    return Region(start, end, context_len, direction=1, cyclic=True)


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


def new_alignment_in_container(container, a, b, type, direction=1):
    query_key = container.alignments[0].query_key

    query_region = container.alignments[0].query_region
    subject_region = container.alignments[0].subject_region
    new_query_region = query_region.new(a, b, allow_wrap=True)
    new_subject_region = subject_region.new(a, b, allow_wrap=True)

    new_subject_region.direction = direction
    alignment = Alignment(new_query_region, new_subject_region, type, query_key, str(uuid4()))
    container.alignments.append(alignment)
    return alignment


test_container = random_container(100, 1, Constants.PCR_PRODUCT)
test_container.alignments = [new_alignment_in_container(test_container, 100, 1000, Constants.PCR_PRODUCT)]

@pytest.fixture(scope='function')
def container():
    return deepcopy(test_container)


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

    groups = AlignmentContainer.redundent_alignment_groups(alignments)
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


class TestExpandPrimers():

    @pytest.mark.parametrize('x', [(75, 90, 1), (2000, 2030, -1)])
    def test_no_expand(self, container, x):
        assert len(container) == 1

        # add primers
        new_alignment_in_container(container, x[0], x[1], Constants.PRIMER, x[2])
        assert len(container) == 2

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 0

    def test_one_pair_two_templates(self, container):
        new_alignment_in_container(container, 110, 900, Constants.PCR_PRODUCT)
        assert len(container) == 2

        # add primers
        new_alignment_in_container(container, 200, 230, Constants.PRIMER)
        new_alignment_in_container(container, 800, 830, Constants.PRIMER, -1)
        assert len(container) == 4

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 6

    def test_one_pair_two_templates2(self, container):
        new_alignment_in_container(container, 1200, 2000, Constants.PCR_PRODUCT)
        assert len(container) == 2

        # add primers
        new_alignment_in_container(container, 200, 230, Constants.PRIMER)
        new_alignment_in_container(container, 1300, 1330, Constants.PRIMER, -1)
        assert len(container) == 4

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 2

    @pytest.mark.parametrize(
        "x", [
            ([(200, 230), (220, 250)], [(750, 800), (800, 830)], 8),
            ([(200, 230), (220, 250)], [(800, 830)], 5),
            ([(200, 230)], [(220, 250), (800, 830)], 5),
            ([(200, 230), (220, 250)], [], 2),
            ([], [(220, 250), (800, 830)], 2)
        ],
        ids=[
            "two_fwd_two_rev",
            "two_fwd_one_rev",
            "one_fwd_two_rev",
            "two_fwd",
            "two_rev"
        ]
    )
    def test_num_alignments_from_pairs(self, x, container):
        # add primers
        for fwd in x[0]:
            new_alignment_in_container(container, fwd[0], fwd[1], Constants.PRIMER)

        for rev in x[1]:
            new_alignment_in_container(container, rev[0], rev[1], Constants.PRIMER, -1)

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == x[2]

    @pytest.mark.parametrize("x", [
        (Constants.PRIMER_MIN_BIND, Constants.PRIMER_MIN_BIND, 3),
        (Constants.PRIMER_MIN_BIND-1, Constants.PRIMER_MIN_BIND, 1),
        (Constants.PRIMER_MIN_BIND, Constants.PRIMER_MIN_BIND-1, 1),
        (Constants.PRIMER_MIN_BIND-1, Constants.PRIMER_MIN_BIND-1, 0),
        ],
        ids=[
            "both primers valid",
            "right primer valid",
            "left primer valid",
            "neither primer valid"
        ]
    )
    def test_expand_one_pair(self, x, container):
        left_primer_len, right_primer_len, expected_num_alignments = x
        assert len(container) == 1

        # add primers
        new_alignment_in_container(container, 100, 100+left_primer_len, Constants.PRIMER)
        new_alignment_in_container(container, 1000-right_primer_len, 1000, Constants.PRIMER, direction=-1)
        assert len(container) == 3

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == expected_num_alignments

    def test_position_of_pairs(self, container):
        left_primer_len, right_primer_len = (Constants.PRIMER_MIN_BIND, Constants.PRIMER_MIN_BIND)
        assert len(container) == 1

        # add primers
        new_alignment_in_container(container, 110, 110+left_primer_len, Constants.PRIMER)
        new_alignment_in_container(container, 900-right_primer_len, 900, Constants.PRIMER, direction=-1)
        assert len(container) == 3

        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 3
        types = set([a.type for a in alignments])
        assert types == set([Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                            Constants.PCR_PRODUCT_WITH_PRIMERS])
        for a in alignments:
            if a.type == Constants.PCR_PRODUCT_WITH_LEFT_PRIMER:
                a.query_region.a == 110
                a.query_region.b == 1000
            elif a.type == Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER:
                a.query_region.a == 100
                a.query_region.b == 900
            elif a.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
                a.query_region.a == 110
                a.query_region.b == 900

    def test_no_expand_if_misses_template(self, container):
        assert len(container) == 1

        # add primers
        new_alignment_in_container(container, 75, 80, Constants.PRIMER)
        new_alignment_in_container(container, 800, 830, Constants.PRIMER, direction=-1)
        assert len(container) == 3

        # expand and check
        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 1

    def test_expand_over_origin(self, container):
        container = random_container(2, 1, Constants.PCR_PRODUCT)
        container.alignments = [new_alignment_in_container(container, 9000, 1000, Constants.PCR_PRODUCT, 1)]

        new_alignment_in_container(container, 9500, 9530, Constants.PRIMER, 1)
        new_alignment_in_container(container, 800, 830, Constants.PRIMER, -1)

        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))

        indices = []
        for a in alignments:
            print(a)
            indices.append((a.query_region.a, a.query_region.b))
        indices.sort()
        print(indices)
        assert indices == [(9000, 830), (9500, 830), (9500, 1000)]

    def test_expand_over_origin2(self, container):
        container = random_container(2, 1, Constants.PCR_PRODUCT)
        container.alignments = [new_alignment_in_container(container, 9000, 1000, Constants.PCR_PRODUCT, 1)]

        new_alignment_in_container(container, 9500, 9530, Constants.PRIMER, -1)
        new_alignment_in_container(container, 800, 830, Constants.PRIMER, 1)

        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 2

        indices = []
        for a in alignments:
            print(a)
            indices.append((a.query_region.a, a.query_region.b))
        indices.sort()
        print(indices)
        assert indices == [(800, 1000), (9000, 9530)]

    def test_overhang_on_primer(self, container):
        new_alignment_in_container(container, 1000-Constants.PRIMER_MIN_BIND, 1100, Constants.PRIMER, -1)
        assert len(container) == 2

        alignments = container.expand_primer_pairs(container.get_groups_by_types(Constants.PCR_PRODUCT))
        assert len(alignments) == 1
        assert alignments[0].query_region.a == 100
        assert alignments[0].query_region.b == 1100


class TestExpandOverlaps():

    def overlap_container(self, container, x1, x2, x3, x4):
        container.alignments = [new_alignment_in_container(container, x1, x2, Constants.PCR_PRODUCT)]
        new_alignment_in_container(container, x3, x4, Constants.PCR_PRODUCT)
        return deepcopy(container)

    @pytest.mark.parametrize('x', [
            (100, 1000, 500, 2000),
            (200, 800, 700, 1000),
            (9000, 1000, 500, 2000)
        ],
        ids=[
            'basic_1',
            'basic_2',
            'over_origin_1']
    )
    def test_expand_overlaps(self, container, x):
        container = self.overlap_container(container, *x)
        assert len(container) == 2

        groups = container.get_groups_by_types(Constants.PCR_PRODUCT)
        alignments = container.expand_overlaps(groups)
        assert len(alignments) == 3

        indices = []
        for a in alignments:
            indices.append((a.query_region.a, a.query_region.b))
        indices.sort()
        assert indices == sorted([(x[0], x[2]), (x[2], x[1]), (x[1], x[3])])

def test_alignments_by_types():
    raise NotImplementedError


def test_types():
    raise NotImplementedError


def test_annotate_fragments():
    raise NotImplementedError
