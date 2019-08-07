from bisect import bisect_left
from typing import List, Dict

import networkx as nx
from Bio.SeqRecord import SeqRecord
from more_itertools import partition, flatten, unique_everseen

from .__version__ import __version__, __title__, __authors__, __homepage__, __repo__
from .blastbiofactory import BioBlastFactory
from .utils import sort_with_keys, bisect_slice_between
from .cost import SpanCost
import itertools
from .utils import Region
import numpy as np
from .log import logger
from shoestring.utils import make_async

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

    PCR_COST = {
        PCR_PRODUCT: 30 + 25 + 10,
        PCR_PRODUCT_WITH_PRIMERS: 30,
        PCR_PRODUCT_WITH_LEFT_PRIMER: 15,
        PCR_PRODUCT_WITH_RIGHT_PRIMER: 15,
        FRAGMENT: 0,
    }

    PRIMER = "PRIMER"  # a primer binding alignment

    MIN_OVERLAP = 20
    MAX_HOMOLOGY = 100
    INF = 10.0 ** 6


class AlignmentException(Exception):
    pass


class Alignment(object):
    """
    A pair of Regions that 'aligns' two regions of DNA sequences.
    """

    __slots__ = ["query_region", "subject_region", "type", "query_key", "subject_key"]

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
        delta_a = qstart - self.query_region.a
        delta_b = qend - self.query_region.b
        subject_copy = self.subject_region.new(
            self.subject_region.a + delta_a,
            self.subject_region.b + delta_b,
            allow_wrap=True,
        )
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


class AlignmentGroup(object):
    """
    A representative Alignment representing a group of alignments.
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


def blast_to_region(query_or_subject, seqdb):
    data = query_or_subject
    record = seqdb[data["origin_key"]]

    s, e = data["start"], data["end"]
    if data["strand"] == -1:
        s, e = e, s
    l = len(record)
    # if data['circular'] and e > l and data['length'] == 2 * l:
    #     e -= l
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
class AlignmentContainer(object):
    valid_types = [
        Constants.PCR_PRODUCT,
        Constants.PRIMER,
        Constants.FRAGMENT,
        Constants.PCR_PRODUCT_WITH_PRIMERS,
        Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
    ]  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord]):
        self.alignments = []
        self.seqdb = seqdb
        self.logger = logger(self)

    def load_blast_json(self, data: dict, type: str):
        """
        Create alignments from a formatted BLAST JSON result.

        :param data: formatted BLAST JSON result
        :param type: the type of alignment to initialize
        :return: None
        """
        self.logger.info("Loading blast json")
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
            self.alignments.append(alignment)

    def annotate_fragments(self, alignments) -> List[Alignment]:
        self.logger.info("annotating fragments")
        annotated = []
        for a in alignments:
            if a.is_perfect_subject() and not a.subject_region.cyclic:
                a.type = Constants.FRAGMENT
                annotated.append(a)
        return annotated

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
        fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.a)
        rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.b)

        pairs = []

        for g in self.logger.tqdm(
            alignment_groups, "INFO", desc="Expanding primer pair"
        ):
            # add products with both existing products
            fwd_bind = bisect_slice_between(
                fwd, fwd_keys, g.query_region.a + 10, g.query_region.b
            )
            rev_bind = bisect_slice_between(
                rev, rev_keys, g.query_region.a, g.query_region.b - 10
            )
            rkeys = [r.query_region.a for r in rev_bind]

            # both primers
            for f in fwd_bind:
                i = bisect_left(rkeys, f.query_region.a)
                for r in rev_bind[i:]:
                    primer_group = g.sub_region(
                        f.query_region.a,
                        r.query_region.b,
                        Constants.PCR_PRODUCT_WITH_PRIMERS,
                    )
                    pairs += primer_group.alignments

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

    def expand_pcr_products(
        self, alignment_groups: List[AlignmentGroup]
    ) -> List[Alignment]:
        """
        Expand the list of alignments from existing regions. Produces new fragments in
        the following two situations:

        ::

            |--------|          alignment 1
                |--------|      alignment 2
            |---|               new alignment

        ::

            |--------|          alignment 1
                 |--------|     alignment 2
                 |---|          new alignment

        :param alignment_groups:
        :return:
        """

        MIN_OVERLAP = Constants.MIN_OVERLAP
        group_sort, group_keys = sort_with_keys(
            alignment_groups, key=lambda x: x.query_region.a
        )
        alignments = []
        for group in logger.tqdm(group_sort, "INFO", desc="expanding pcr products"):
            i = bisect_left(group_keys, group.query_region.a)
            arr, keys = sort_with_keys(group_sort[i:], key=lambda x: x.query_region.a)

            # get list of overlapping alignments
            overlapping = bisect_slice_between(
                arr, keys, group.query_region.a, group.query_region.b
            )

            #
            for og in overlapping:
                if og is not group:
                    if og.query_region.a - group.query_region.a > MIN_OVERLAP:
                        ag1 = group.sub_region(
                            group.query_region.a, og.query_region.a, Constants.PCR_PRODUCT
                        )
                        alignments += ag1.alignments
                    if group.query_region.b - og.query_region.a > MIN_OVERLAP:
                        ag2 = group.sub_region(
                            og.query_region.a, group.query_region.b, Constants.PCR_PRODUCT
                        )
                        alignments += ag2.alignments
        return alignments

    # TODO: expand should just add more
    def expand(self, expand_overlaps=True, expand_primers=True):

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
            expanded = self.expand_pcr_products(templates)
            self.alignments += expanded
            self.logger.info("Number of new alignments: {}".format(len(expanded)))
        self.logger.info("Number of total alignments: {}".format(len(self.alignments)))
        self.logger.info(
            "Number of total groups: {}".format(len(self.alignment_groups))
        )

    # TODO: change 'start' and 'end' to left and right end for regions...

    # TODO: type increases computation time exponentially, we only need 'b' from group1 and 'a' from group2
    @staticmethod
    def alignment_hash(a):
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
            grouped.setdefault(cls.alignment_hash(a), list()).append(a)
        alignment_groups = []
        for group in grouped.values():
            alignment_groups.append(AlignmentGroup(group[0].query_region, group))
        return alignment_groups

    @property
    def types(self) -> List[str]:
        """
        Return all valid types

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


class AssemblyGraphBuilder(object):
    def __init__(self, alignment_container: AlignmentContainer, span_cost=None):
        self.container = alignment_container
        if span_cost is None:
            self.span_cost = SpanCost()
        else:
            self.span_cost = span_cost
        self.G = None
        self.logger = logger(self)

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
            q = g.query_region
            b_expand = True
            a_expand = True

            # TODO: Constants should be a class?
            if g.type in [Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
                b_expand = False
            if g.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
                a_expand = False

            ### INTERNAL EDGE
            if g.type == Constants.FRAGMENT:
                internal_cost = 0
            elif g.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER]:
                internal_cost = 30 + 30
            elif g.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
                internal_cost = 30
            elif g.type == Constants.PCR_PRODUCT:
                internal_cost = 30 + 30 + 30

            self.G.add_edge(
                (q.a, a_expand, 'A'),
                (q.b, b_expand, 'B'),
                weight=internal_cost,
                name='',
                span_length=len(g.query_region)
            )

            a_arr.add((q.a, a_expand, 'A'))
            b_arr.add((q.b, b_expand, 'B'))

        ### EXTERNAL EDGES
        if groups:
            query = groups[0].query_region
            for (b, b_expand, bid), (a, a_expand, aid) in itertools.product(b_arr, a_arr):
                if a != b:
                    ba = query.new(b, a)
                    ab = query.new(a, b)

                    # TODO: no way to determine overlaps from just end points

                    r = ba # sorted([(r1, len(r1)), (r2, len(r2))], key=lambda x: x[1])[0][0]
                    cost = self.span_cost.cost(len(r), (b_expand, a_expand))
                    if cost < 10000:
                        self.G.add_edge(
                            (b, b_expand, bid),
                            (a, a_expand, aid),
                            weight=cost,
                            name='',
                            span_length=len(r)
                        )
        else:
            self.logger.warn("There is nothing to assembly. There are no alignments.")
        return self.G
