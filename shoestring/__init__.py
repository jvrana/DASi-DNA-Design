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


# TODO: having the costs mixed up between classes is confusing
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
    ]

    def __init__(self, seqdb: Dict[str, SeqRecord]):
        self.alignments = []
        self.seqdb = seqdb

    def load_blast_json(self, data: dict, type: str):
        """
        Create alignments from a formatted BLAST JSON result.

        :param data: formatted BLAST JSON result
        :param type: the type of alignment to initialize
        :return: None
        """
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
        annotated = []
        for a in alignments:
            if a.is_perfect_subject() and not a.subject_region.circular:
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

        for g in alignment_groups:
            # add products with both existing products
            fwd_bind = bisect_slice_between(
                fwd, fwd_keys, g.query_region.a + 10, g.query_region.b
            )
            rev_bind = bisect_slice_between(
                rev, rev_keys, g.query_region.a, g.query_region.b - 10
            )
            rkeys = [r.query_region.a for r in rev_bind]
            for f in fwd_bind:
                i = bisect_left(rkeys, f.query_region.a)
                for r in rev_bind[i:]:
                    primer_group = g.sub_region(
                        f.query_region.a,
                        r.query_region.b,
                        Constants.PCR_PRODUCT_WITH_PRIMERS,
                    )
                    pairs += primer_group.alignments

            # add products with one existing primer
            for f in fwd_bind:
                left_primer_group = g.sub_region(
                    f.query_region.a,
                    g.query_region.b,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                )
                pairs += left_primer_group.alignments
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
        group_sort, group_keys = sort_with_keys(
            alignment_groups, key=lambda x: x.query_region.a
        )
        alignments = []
        for g in group_sort:
            i = bisect_left(group_keys, g.query_region.a)
            arr, keys = sort_with_keys(group_sort[i:], key=lambda x: x.query_region.a)
            overlapping = bisect_slice_between(
                arr, keys, g.query_region.a, g.query_region.b
            )

            for og in overlapping:
                if og is not g:
                    if og.query_region.a - g.query_region.a > 20:
                        ag1 = g.sub_region(
                            g.query_region.a, og.query_region.a, Constants.PCR_PRODUCT
                        )
                        alignments += ag1.alignments
                    if g.query_region.b - og.query_region.a > 20:
                        ag2 = g.sub_region(
                            og.query_region.a, g.query_region.b, Constants.PCR_PRODUCT
                        )
                        alignments += ag2.alignments
        return alignments

    def expand(self):
        print()
        print("=== Expanding alignments ===")
        # We annotate any original PCR_PRODUCT with FRAGMENT if they are 'perfect_subjects'
        # This means they already exist as pre-made fragments
        print("Number of alignments: {}".format(len(self.alignments)))
        annotated = self.annotate_fragments(
            self.get_alignments_by_types(Constants.PCR_PRODUCT)
        )

        print("Number of perfect subjects: {}".format(len(annotated)))
        templates = self.get_groups_by_types(
            [Constants.PCR_PRODUCT, Constants.FRAGMENT]
        )

        pairs = self.expand_primer_pairs(templates)
        self.alignments += pairs
        print("Number of pairs: {}".format(len(pairs)))

        expanded = self.expand_pcr_products(templates)
        self.alignments += expanded
        print("Number of new alignments: {}".format(len(expanded)))
        print("Number of total alignments: {}".format(len(self.alignments)))
        print("Number of total groups: {}".format(len(self.alignment_groups)))

    # TODO: change 'start' and 'end' to left and right end for regions...

    @staticmethod
    def alignment_hash(a):
        return (a.query_region.a, a.query_region.b, a.query_region.direction, a.type)

    @classmethod
    def group(cls, alignments):
        grouped = {}
        for a in alignments:
            grouped.setdefault(cls.alignment_hash(a), list()).append(a)
        alignment_groups = []
        for group in grouped.values():
            alignment_groups.append(AlignmentGroup(group[0].query_region, group))
        return alignment_groups

    @property
    def types(self):
        return tuple(self.valid_types)

    def get_groups_by_types(self, types: List[str]) -> List[AlignmentGroup]:
        groups = self.groups_by_type
        if isinstance(types, str):
            return groups[types]
        else:
            return list(unique_everseen(flatten([groups[t] for t in types])))

    def get_alignments_by_types(self, types: List[str]) -> List[Alignment]:
        groups = self.get_groups_by_types(types)
        return list(flatten([g.alignments for g in groups]))

    @property
    def groups_by_type(self):
        d = {}
        for t in self.valid_types:
            d[t] = []
        for a in self.alignments:
            d[a.type].append(a)
        return {t: self.group(v) for t, v in d.items()}

    @property
    def alignment_groups(self):
        return self.group(self.alignments)


class AssemblyGraphBuilder(object):
    def __init__(self, alignment_container: AlignmentContainer):
        self.container = alignment_container
        self.span_cost = SpanCost()
        self.G = None

    def _edge_weight(self, type1: str, type2: str, span: int) -> float:
        ext = self._extension_tuple_from_types(type1, type2)
        assembly_cost = np.clip(self.span_cost.cost(span, ext), 0, Constants.INF)
        frag_cost1 = Constants.PCR_COST[type1] # we only count the source fragment cost
        return assembly_cost + frag_cost1

    def _extension_tuple_from_types(self, left_type, right_type):
        """
        Converts two "TYPES" into an 'extension tuple' that can be consumed by
        the cost classes. The 'extension tuple' is a tuple representing the ability
        of the fragments across a junction to be design to bridge a gap. For example:

        (0, 0) means the following gap cannot be spanned:

        ::
            (0, 0)
            ------|      |------

        It is also important to note that "0" means the either the primer for that fragment
        exists

        While the follow gap *may* be spanned by designing the left fragment with a longer
        primer:

        ::

            (1, 0)
            -------|    |-------
               <----------

        And finally, the following gap may be spanned by both primers.


        ::

            (1, 1)
                    -------->
            -------|    |-------
               <----------

        :param left_type:
        :param right_type:
        :return:
        """
        if (
            left_type == Constants.PCR_PRODUCT_WITH_PRIMERS
            or Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER
        ):
            left_extendable = 1
        else:
            left_extendable = 0

        if (
            right_type == Constants.PCR_PRODUCT_WITH_PRIMERS
            or Constants.PCR_PRODUCT_WITH_LEFT_PRIMER
        ):
            right_extendable = 1
        else:
            right_extendable = 0
        return (left_extendable, right_extendable)

    def _add_node(self, g1: AlignmentGroup):
        """Add a node from an alignment group to the graph."""
        n1 = self.container.alignment_hash(g1.alignments[0])
        self.G.add_node(
            n1,
            start=g1.query_region.start,
            end=g1.query_region.end,
            region=str(g1.query_region),
        )
        return n1

    def _add_edge(self, g1, g2, span_length, **kwargs):
        """Add nodes, edge, and edge_weight from two alignment groups. A spanning distance (bp)
        must be provided."""
        n1 = self._add_node(g1)
        n2 = self._add_node(g2)
        edata = {
            "weight": self._edge_weight(g1.type, g2.type, span_length),
            "region_1": str(g1.query_region),
            "region_2": str(g2.query_region),
        }
        edata.update(kwargs)
        self.G.add_edge(n1, n2, **edata)

    def add_edge(self, group, other_group):
        """
        Create an edge between one alignment group and another.

        :param group:
        :param other_group:
        :return:
        """
        # verify contexts are the same
        r1 = group.query_region
        r2 = other_group.query_region
        assert r2.same_context(r1)

        if r1 in r2 or r2 in r1:
            return
        if r2.a in r1 and r2.b not in r1:
            # strict forward overlap
            overlap = r1.intersection(r2)
            if not overlap:
                raise Exception("We expected an overlap here")
            # TODO: penalize small overlap
            self._add_edge(
                group, other_group, span_length=-len(overlap), name="overlap"
            )
        elif r2.a not in r1 and r2.b not in r1:
            try:
                connecting_span = r1.connecting_span(r2)
            except Exception as e:
                raise e
            if connecting_span:
                span_length = len(connecting_span)
            elif r1.consecutive(r2):
                span_length = 0
            else:
                return
            self._add_edge(
                group, other_group, span_length=span_length, name="synthesis"
            )

    # TODO: linear assemblies
    # TODO: replace region methods with something faster.
    # TODO: edge cost should include (i) cost of producing the source and (ii) cost of assembly
    # Verify queries have same context
    def build_assembly_graph(self):
        self.G = nx.DiGraph(name="Assembly Graph")
        # TODO: tests for validating over-origin edges are being produced
        # produce non-spanning edges
        # makes edges for any regions
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
        for g1, g2 in itertools.product(groups, repeat=2):
            self.add_edge(g1, g2)
        return self.G
