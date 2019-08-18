from bisect import bisect_left
from typing import List, Dict

import networkx as nx
from Bio.SeqRecord import SeqRecord
from more_itertools import partition, flatten, unique_everseen

from .__version__ import __version__, __title__, __authors__, __homepage__, __repo__
from .utils import sort_with_keys, bisect_slice_between
from .cost import SpanCost
import itertools
from .utils import Region
from .log import logger
from dasi.utils import make_async
from collections.abc import Sized


# TODO: instead of matching fragments, match 'b's to 'a's (e.g. 360000 edges vs 6390 edges!)
# TODO: having the costs mixed up between classes is confusing
# TODO: move these to their own classes?

class Constants(object):
    FRAGMENT = (
        "PRE-MADE DNA FRAGMENT"
    )  # an alignment that is generate from an already existing PCR product or fragment
    PCR_PRODUCT = (
        "PCR_PRODUCT"
    )  # an alignment that is to be generated from a PCR product
    PCR_PRODUCT_WITH_PRIMERS = (
        "PCR_PRODUCT_WITH_PRIMERS"
    )  # PCR product that can be produces from existing primers
    PCR_PRODUCT_WITH_LEFT_PRIMER = (
        "PCR_PRODUCT_WITH_LEFT_PRIMER"
    )  # PCR product with existing left primer
    PCR_PRODUCT_WITH_RIGHT_PRIMER = (
        "PCR_PRODUCT_WITH_RIGHT_PRIMER"
    )  # PCR product with existing right primer
    SHARED_FRAGMENT = (
        "FRAGMENT_SHARED_WITH_OTHER_QUERIES"
    )  # A fragment alignment that is shared with other queries for potential reuse

    PCR_COST = {
        PCR_PRODUCT: 30 + 25 + 10,
        PCR_PRODUCT_WITH_PRIMERS: 30,
        PCR_PRODUCT_WITH_LEFT_PRIMER: 15,
        PCR_PRODUCT_WITH_RIGHT_PRIMER: 15,
        FRAGMENT: 0,
    }

    PRIMER = "PRIMER"  # a primer binding alignment

    PRIMER_MIN_BIND = 15
    MIN_OVERLAP = 20
    MAX_HOMOLOGY = 100
    INF = 10.0 ** 6


class DASiException(Exception):
    pass


class AlignmentException(DASiException):
    pass


class AlignmentContainerException(DASiException):
    pass


ALIGNMENT_SLOTS = ["query_region", "subject_region", "type", "query_key", "subject_key"]


class Alignment(object):
    """
    A pair of Regions that 'aligns' two regions of DNA sequences. Regions must be
    the same length.

    A subregion of both regions may be taken.
    """

    __slots__ = ALIGNMENT_SLOTS[:]

    def __init__(
        self,
        query_region: Region,
        subject_region: Region,
        type: str,
        query_key: str,
        subject_key: str,
    ):
        """
        Makes an alignment between two regions of sequences. Validates the regions are the same length.

        :param query_region: Query region this alignment aligns to
        :param subject_region: Subject region this alignment aligns to.
        :param type: Type of alignment
        :param query_key: The record identifier for the query
        :param subject_key: The record identifier for the subject
        """
        self.query_region = query_region
        self.subject_region = subject_region
        self.validate()
        self.type = type
        self.query_key = query_key
        self.subject_key = subject_key

    def validate(self):
        if not len(self.query_region) == len(self.subject_region):
            raise AlignmentException(
                "Regions must have the same size: {} vs {}".format(
                    len(self.query_region), len(self.subject_region)
                )
            )

    def is_perfect_subject(self):
        return len(self.subject_region) == self.subject_region.context_length

    def sub_region(self, qstart: int, qend: int, type=None):
        """
        Returns a copy of the alignment between the inclusive start and end relative to the
        query region.

        :param qstart: start of the query sub region
        :param qend: end of the query sub region
        :param type: optional type of alignment to return
        :return:
        """
        query_copy = self.query_region.sub(qstart, qend)
        if self.query_region.a != qstart:
            delta_a_span = self.query_region.sub(self.query_region.a, qstart)
            delta_a = len(delta_a_span)
        else:
            delta_a = 0
        if qend != self.query_region.b:
            delta_b_span = self.query_region.sub(qend, self.query_region.b)
            delta_b = -len(delta_b_span)
        else:
            delta_b = 0
        if delta_a == 0:
            delta_a = None
        if delta_b == 0:
            delta_b = None
        subject_copy = self.subject_region[delta_a:delta_b]

        if type is None:
            type = self.type
        self.validate()
        return self.__class__(
            query_region=query_copy,
            subject_region=subject_copy,
            type=type,
            query_key=self.query_key,
            subject_key=self.subject_key,
        )

    def __str__(self):
        return "<{} {} {} {}>".format(
            self.__class__.__name__, self.type, self.query_region, self.subject_region
        )

    def __repr__(self):
        return str(self)

class AlignmentGroup(object):
    """
    A representative Alignment representing a group of alignments sharing the
    same starting and ending position for a query sequence.
    """

    __slots__ = ["query_region", "alignments", "name"]

    def __init__(self, query_region: Region, alignments: List[Alignment], name=None):
        self.query_region = query_region
        self.alignments = alignments
        self.name = name

    @property
    def type(self):
        return self.alignments[0].type

    # TODO: making subregions for all alignments takes a LONG TIME (5X shorter if you skip this).
    def sub_region(self, qstart: int, qend: int, type: str):
        alignments_copy = [a.sub_region(qstart, qend) for a in self.alignments]
        for a in alignments_copy:
            a.type = type
        return self.__class__(
            query_region=self.query_region.sub(qstart, qend),
            alignments=alignments_copy,
            name="subregion",
        )


# TODO: make a test for this conversion
def blast_to_region(query_or_subject, seqdb):
    data = query_or_subject
    record = seqdb[data["origin_key"]]

    s, e = data["start"], data["end"]
    l = len(record)
    # if data['circular'] and e > l and data['length'] == 2 * l:
    #     e -= l
    if data['strand'] == -1:
        s, e = e, s
    region = Region(
        s - 1,
        e,
        length=l,
        cyclic=data["circular"],
        direction=data["strand"],
        index=0,
        name="{}: {}".format(record.id, record.name),
        allow_wrap=True,
    )
    return region


# TODO: This assumes a single query. Verify this.
# TODO: need way to run blast from multiple databases on a single query
class AlignmentContainer(Sized):
    """
    Container for a set of query-to-subject alignments for a single query.

    This class contains:
        1. Methods grouping alignments together according to those
    alignments that share the same starting and ending points.
        2. Methods for 'expanding' alignments
    """

    valid_types = (
        Constants.PCR_PRODUCT,
        Constants.PRIMER,
        Constants.FRAGMENT,
        Constants.PCR_PRODUCT_WITH_PRIMERS,
        Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
        Constants.SHARED_FRAGMENT
    )  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord], alignments=None):
        self._alignments = []
        if alignments is not None:
            self.alignments = alignments
        self.seqdb = seqdb
        self.logger = logger(self)

    @property
    def alignments(self):
        return self._alignments

    @alignments.setter
    def alignments(self, v):
        self._alignments = v
        self._check_single_query_key(self._alignments)

    @staticmethod
    def _check_single_query_key(alignments):
        keys = set(a.query_key for a in alignments)
        if len(keys) > 1:
            raise AlignmentContainerException("AlignmentContainer cannot contain more than one query. Contains the following"
                             "query keys: {}".format(keys))

    def annotate_fragments(self, alignments) -> List[Alignment]:
        self.logger.info("annotating fragments")
        annotated = []
        for a in alignments:
            if a.is_perfect_subject() and not a.subject_region.cyclic:
                a.type = Constants.FRAGMENT
                annotated.append(a)
        return annotated

    # TODO: test for this method
    @classmethod
    def filter_alignments_by_span(cls, alignments, region, key=None):
        fwd, fwd_keys = sort_with_keys(alignments, key=key)
        found = []
        for a, b in region.ranges():
            found += bisect_slice_between(
                fwd, fwd_keys, a, b
            )
        return found

    def expand_primer_pairs(
        self, alignment_groups: List[AlignmentGroup]
    ) -> List[Alignment]:
        """
        Creates new alignments for all possible primer pairs. Searches for fwd and
        rev primer pairs that exist within other alignments and produces all combinations
        of alignments that can form from these primer pairs.

        :return: list
        """
        primers = self.get_alignments_by_types(Constants.PRIMER)

        rev, fwd = partition(lambda p: p.subject_region.direction == 1, primers)
        # fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.b)
        # rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.a)

        pairs = []

        for g in self.logger.tqdm(
            alignment_groups, "INFO", desc="Expanding primer pair"
        ):
            query_ranges = g.query_region.ranges()
            fwd_bind_region = g.query_region[Constants.PRIMER_MIN_BIND:]
            rev_bind_region = g.query_region[:-Constants.PRIMER_MIN_BIND]
            fwd_bind = self.filter_alignments_by_span(fwd, fwd_bind_region, key=lambda p: p.query_region.b)
            rev_bind = self.filter_alignments_by_span(rev, rev_bind_region, key=lambda p: p.query_region.a)
            rkeys = [r.query_region.a for r in rev_bind]

            # both primers
            for f in fwd_bind:
                _rev_bind = []
                if len(query_ranges) == 1:
                    i = bisect_left(rkeys, f.query_region.a)
                    _rev_bind = rev_bind[i:]
                else:
                    _rev_span = g.query_region.sub(f.query_region.a, g.query_region.b)
                    for a, b in _rev_span.ranges():

                        _rev_bind += bisect_slice_between(
                            rev_bind, rkeys, a, b
                        )
                for r in _rev_bind:
                    try:
                        primer_group = g.sub_region(
                            f.query_region.a,
                            r.query_region.b,
                            Constants.PCR_PRODUCT_WITH_PRIMERS,
                        )
                        pairs += primer_group.alignments
                    except IndexError:
                        pass

            # left primer
            for f in fwd_bind:
                left_primer_group = g.sub_region(
                    f.query_region.a,
                    g.query_region.b,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                )
                pairs += left_primer_group.alignments

            # right primer
            for r in rev_bind:
                right_primer_group = g.sub_region(
                    g.query_region.a,
                    r.query_region.b,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                )
                pairs += right_primer_group.alignments
        return pairs

    def expand_overlaps(
        self, alignment_groups: List[AlignmentGroup],
            type=Constants.PCR_PRODUCT
    ) -> List[Alignment]:
        """
        Expand the list of alignments from existing regions. Produces new fragments in
        the following three situations:

        ::

            |--------|          alignment 1
                |--------|      alignment 2
            |---|               new alignment


            |--------|          alignment 1
                 |--------|     alignment 2
                 |---|          new alignment


            |--------|          alignment 1
                 |--------|     alignment 2
                     |----|     new alignment

        :param alignment_groups:
        :return: list
        """

        MIN_OVERLAP = Constants.MIN_OVERLAP
        group_sort, group_keys = sort_with_keys(
            alignment_groups, key=lambda x: x.query_region.a
        )
        alignments = []
        for group_a in logger.tqdm(group_sort, "INFO", desc="expanding pcr products"):
            i = bisect_left(group_keys, group_a.query_region.a)
            arr, keys = sort_with_keys(group_sort[i:], key=lambda x: x.query_region.a)

            # get list of overlapping alignments
            overlapping = bisect_slice_between(
                arr, keys, group_a.query_region.a, group_a.query_region.b
            )

            #
            for group_b in overlapping:
                if group_b is not group_a:
                    if group_b.query_region.a - group_a.query_region.a > MIN_OVERLAP:
                        ag1 = group_a.sub_region(
                            group_a.query_region.a, group_b.query_region.a, type
                        )
                        alignments += ag1.alignments
                    if group_a.query_region.b - group_b.query_region.a > MIN_OVERLAP:
                        ag2 = group_a.sub_region(
                            group_b.query_region.a, group_a.query_region.b, type
                        )
                        alignments += ag2.alignments
                    if group_b.query_region.b - group_a.query_region.b > MIN_OVERLAP:
                        ag3 = group_b.sub_region(
                            group_a.query_region.b, group_b.query_region.b, type
                        )
                        alignments += ag3.alignments
        return alignments

    # TODO: expand should just add more
    def expand(self, expand_overlaps=True, expand_primers=True):
        """
        Expand the number of alignments in this container using overlaps or primers.

        :param expand_overlaps:
        :param expand_primers:
        :return:
        """

        self.logger.info("=== Expanding alignments ===")
        # We annotate any original PCR_PRODUCT with FRAGMENT if they are 'perfect_subjects'
        # This means they already exist as pre-made fragments
        self.logger.info("Number of alignments: {}".format(len(self.alignments)))

        # TODO: what is annotate_fragments for?
        annotated = self.annotate_fragments(
            self.get_alignments_by_types(Constants.PCR_PRODUCT)
        )

        self.logger.info("Number of perfect subjects: {}".format(len(annotated)))
        templates = self.get_groups_by_types(
            [Constants.PCR_PRODUCT, Constants.FRAGMENT]
        )

        if expand_primers:
            pairs = self.expand_primer_pairs(templates)
            self.alignments += pairs
            self.logger.info("Number of pairs: {}".format(len(pairs)))

        if expand_overlaps:
            expanded = self.expand_overlaps(templates)
            self.alignments += expanded
            self.logger.info("Number of new alignments: {}".format(len(expanded)))
        self.logger.info("Number of total alignments: {}".format(len(self.alignments)))
        self.logger.info(
            "Number of total groups: {}".format(len(self.alignment_groups))
        )

    # TODO: change 'start' and 'end' to left and right end for regions...
    # TODO: type increases computation time exponentially, we only need 'b' from group1 and 'a' from group2
    @staticmethod
    def _alignment_hash(a):
        """A hashable representation of an alignment for grouping."""
        return (a.query_region.a, a.query_region.b, a.query_region.direction, a.type)

    @classmethod
    def group(cls, alignments: List[Alignment]) -> List[AlignmentGroup]:
        """
        Return AlignmentGroups that have been grouped by alignment_hash

        :param alignments:
        :return:
        """
        grouped = {}
        for a in alignments:
            grouped.setdefault(cls._alignment_hash(a), list()).append(a)
        alignment_groups = []
        for group in grouped.values():
            alignment_groups.append(AlignmentGroup(group[0].query_region, group))
        return alignment_groups

    @property
    def types(self) -> List[str]:
        """
        Return all valid types.

        :return:
        """
        return tuple(self.valid_types)

    def get_groups_by_types(self, types: List[str]) -> List[AlignmentGroup]:
        """
        Return AlignmentGroups by fragment type

        :param types: list of types
        :return:
        """
        groups = self.groups_by_type
        if isinstance(types, str):
            return groups[types]
        else:
            return list(unique_everseen(flatten([groups[t] for t in types])))

    def get_alignments_by_types(self, types: List[str]) -> List[Alignment]:
        groups = self.get_groups_by_types(types)
        return list(flatten([g.alignments for g in groups]))

    @property
    def groups_by_type(self) -> Dict[str, AlignmentGroup]:
        """
        Return alignment groups according to their alignment 'type'

        :return: dict
        """
        d = {}
        for t in self.valid_types:
            d[t] = []
        for a in self.alignments:
            d[a.type].append(a)
        return {t: self.group(v) for t, v in d.items()}

    @property
    def alignment_groups(self):
        """
        Return all alignment groups

        :return:
        """
        return self.group(self.alignments)

    def __len__(self):
        return len(self.alignments)


class AlignmentContainerFactory(object):
    """
    Class that maintains a shared list of alignments and shared sequence database.

    AlignmentContainers can be retrieved in a dict grouped by their query via `.containers()`
    """
    valid_types = (
        Constants.PCR_PRODUCT,
        Constants.PRIMER,
        Constants.FRAGMENT,
        Constants.PCR_PRODUCT_WITH_PRIMERS,
        Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
        Constants.SHARED_FRAGMENT
    )  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord]):
        self.alignments = {}
        self.logger = logger(self)
        self.seqdb = seqdb

    def load_blast_json(self, data: dict, type: str):
        """
        Create alignments from a formatted BLAST JSON result.

        :param data: formatted BLAST JSON result
        :param type: the type of alignment to initialize
        :return: None
        """
        self.logger.info("Loading blast json ({} entries) to fragment type \"{}\"".format(len(data), type))
        assert type in self.valid_types
        for d in data:
            query_region = blast_to_region(d["query"], self.seqdb)
            subject_region = blast_to_region(d["subject"], self.seqdb)
            query_key = d["query"]["origin_key"]
            subject_key = d["subject"]["origin_key"]

            alignment = Alignment(
                query_region,
                subject_region,
                type=type,
                query_key=query_key,
                subject_key=subject_key,
            )
            self.alignments.setdefault(query_key, list()).append(alignment)

    def containers(self):
        container_dict = {}
        for key, alignments in self.alignments.items():
            container_dict[key] = AlignmentContainer(self.seqdb, alignments=alignments)
        return container_dict


class AssemblyGraphBuilder(object):
    """
    Class that builds an AssemblyGraph from an alignment container.
    """

    COST_THRESHOLD = 10000

    def __init__(self, alignment_container: AlignmentContainer, span_cost=None):
        self.container = alignment_container
        if span_cost is None:
            self.span_cost = SpanCost()
        else:
            self.span_cost = span_cost
        self.G = None
        self.logger = logger(self)

    def alignment_to_internal_edge(self, align):
        q = align.query_region
        a_expand, b_expand = True, True
        if align.type in [Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            b_expand = False
        if align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            a_expand = False

        ### INTERNAL EDGE
        if align.type == Constants.FRAGMENT:
            internal_cost = 0
        elif align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER]:
            internal_cost = 30 + 30
        elif align.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
            internal_cost = 30
        elif align.type == Constants.PCR_PRODUCT:
            internal_cost = 30 + 30 + 30

        return (q.a, a_expand, 'A'), (q.b, b_expand, 'B'), dict(
                    weight=internal_cost,
                    name='',
                    span_length=len(align.query_region),
                    type=align.type
                )

    def build_assembly_graph(self):

        self.G = nx.DiGraph(name="Assembly Graph")

        groups, group_keys = sort_with_keys(
            self.container.get_groups_by_types(
                [
                    Constants.PCR_PRODUCT,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_PRIMERS,
                    Constants.FRAGMENT,
                ]
            ),
            key=lambda g: g.query_region.a,
        )

        # TODO: reduce number of times edge

        a_arr = set()
        b_arr = set()

        for g in groups:
            a, b, ab_data = self.alignment_to_internal_edge(g)
            self.G.add_edge(a, b, **ab_data)
            a_arr.add(a)
            b_arr.add(b)

        ### EXTERNAL EDGES
        if groups:
            query = groups[0].query_region
            for (b, b_expand, bid), (a, a_expand, aid) in itertools.product(b_arr, a_arr):
                if a != b:
                    ba = query.new(b, a)
                    ab = query.new(a, b)

                    # TODO: PRIORITY no way to determine overlaps from just end points

                    r = ba # sorted([(r1, len(r1)), (r2, len(r2))], key=lambda x: x[1])[0][0]
                    cost, desc = self.span_cost.cost_and_desc(len(r), (b_expand, a_expand))
                    if cost < self.COST_THRESHOLD:
                        self.G.add_edge(
                            (b, b_expand, bid),
                            (a, a_expand, aid),
                            weight=cost,
                            name='',
                            span_length=len(r),
                            type=desc
                        )
        else:
            self.logger.warn("There is nothing to assembly. There are no alignments.")
        return self.G
