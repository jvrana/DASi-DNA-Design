from .__version__ import __version__, __title__, __authors__, __homepage__, __repo__

from lobio.models import Region, Context
import networkx as nx
from .utils import sort_with_keys, bisect_slice_between
from bisect import bisect_left
from typing import List, Dict
from Bio.SeqRecord import SeqRecord
from more_itertools import partition, flatten, unique_everseen
from .blastbiofactory import BioBlastFactory


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
    PRIMER = "PRIMER"  # a primer binding alignment

    COLOR = "edge_type"
    RED = "molecule"
    BLUE = "assembly"

    MAX_HOMOLOGY = 100


class AlignmentException(Exception):
    pass


class Alignment(object):

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
        return len(self.subject_region) == len(self.subject_region.context)

    def sub_region(self, qstart: int, qend: int, type=None):
        """
        Returns a copy of the alignment between the inclusive start and end relative to the
        query region.

        :param qstart: start of the query sub region
        :param qend: end of the query sub region
        :param type: optional type of alignment to return
        :return:
        """
        query_copy = self.query_region.sub_region(qstart, qend)
        subject_copy = self.subject_region.copy()
        delta_s = qstart - self.query_region.left_end
        delta_e = self.query_region.right_end - qend
        subject_copy.extend_left_end(-delta_s)
        subject_copy.extend_right_end(-delta_e)
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

    __slots__ = ["query_region", "alignments"]

    def __init__(self, query_region: Region, alignments: List[Alignment]):
        self.query_region = query_region
        self.alignments = alignments

    def sub_region(self, qstart: int, qend: int, type: str):
        return self.__class__(
            query_region=self.query_region.sub_region(qstart, qend),
            alignments=[a.sub_region(qstart, qend, type) for a in self.alignments],
        )


def blast_to_region(query_or_subject, seqdb):
    data = query_or_subject
    l = data["length"]
    record = seqdb[data["origin_key"]]
    l = len(record)
    # if data['circular']:
    #     assert l % 2 == 0
    #     l = int(l / 2.0)
    s, e = data["start"], data["end"]
    if data["strand"] == -1:
        s, e = e, s
    region = Region(
        s,
        e,
        direction=data["strand"],
        context=Context(l, data["circular"], start_index=1),
        name="{}: {}".format(record.id, record.name),
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
        fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.left_end)
        rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.right_end)

        pairs = []

        for g in alignment_groups:
            # add products with both existing products
            fwd_bind = bisect_slice_between(
                fwd, fwd_keys, g.query_region.left_end + 10, g.query_region.right_end
            )
            rev_bind = bisect_slice_between(
                rev, rev_keys, g.query_region.left_end, g.query_region.right_end - 10
            )
            rkeys = [r.query_region.left_end for r in rev_bind]
            for f in fwd_bind:
                i = bisect_left(rkeys, f.query_region.left_end)
                for r in rev_bind[i:]:
                    primer_group = g.sub_region(
                        f.query_region.left_end,
                        r.query_region.right_end,
                        Constants.PCR_PRODUCT_WITH_PRIMERS,
                    )
                    pairs += primer_group.alignments

            # add products with one existing primer
            for f in fwd_bind:
                left_primer_group = g.sub_region(
                    f.query_region.left_end,
                    g.query_region.right_end,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                )
                pairs += left_primer_group.alignments
            for r in rev_bind:
                right_primer_group = g.sub_region(
                    g.query_region.left_end,
                    r.query_region.right_end,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                )
                pairs += right_primer_group.alignments
        return pairs

    def expand_pcr_products(
        self, alignment_groups: List[AlignmentGroup]
    ) -> List[Alignment]:
        group_sort, group_keys = sort_with_keys(
            alignment_groups, key=lambda x: x.query_region.left_end
        )
        alignments = []
        for g in group_sort:
            i = bisect_left(group_keys, g.query_region.left_end)
            arr, keys = sort_with_keys(
                group_sort[i:], key=lambda x: x.query_region.left_end
            )
            overlapping = bisect_slice_between(
                arr, keys, g.query_region.left_end, g.query_region.right_end
            )

            for og in overlapping:
                if og is not g:
                    if og.query_region.left_end - g.query_region.left_end > 20:
                        ag1 = g.sub_region(
                            g.query_region.left_end,
                            og.query_region.left_end,
                            Constants.PCR_PRODUCT,
                        )
                        alignments += ag1.alignments
                    if g.query_region.right_end - og.query_region.left_end > 20:
                        ag2 = g.sub_region(
                            og.query_region.left_end,
                            g.query_region.right_end,
                            Constants.PCR_PRODUCT,
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

    # Verify queries have same context
    def build_assembly_graph(self):
        G = nx.DiGraph(name="Assembly Graph")
        self.expand()

        def add_edge(start, stop, length, color):
            G.add_edge(
                "{}_{}".format(start, color),
                "{}_{}".format(stop, color),
                **{Constants.COLOR: color, "length": length, "weight": 1}
            )

        # RED edges
        for g in self.alignment_groups:
            add_edge(
                g.query_region.left_end,
                g.query_region.right_end,
                length=len(g.query_region),
                color=Constants.RED,
            )

        # BLUE edges
        # TODO: tests for validating over-origin edges are being produced
        groups, group_keys = sort_with_keys(
            self.alignment_groups, key=lambda g: g.query_region.left_end
        )
        for g in groups:

            try:
                homology = g.query_region[-Constants.MAX_HOMOLOGY :]

                i = bisect_left(group_keys, homology.left_end)
                other_groups = groups[i:]
                if homology.spans_origin():
                    i = bisect_left(group_keys, homology.right_end)
                    other_groups += groups[:i]

                for g2 in other_groups:
                    if g2 is not g:
                        if not g2.query_region.context == g.query_region.context:
                            assert g2.query_region.context == g.query_region.context

                        if g.query_region.encompasses(g2.query_region):
                            continue
                        else:
                            overlap = g.query_region.get_overlap(g2.query_region)
                            if not overlap:
                                overlap = Region(
                                    g.query_region.right_end,
                                    g2.query_region.left_end,
                                    context=g.query_region.context,
                                )
                            add_edge(
                                g.query_region.right_end,
                                g2.query_region.left_end,
                                length=len(overlap),
                                color=Constants.BLUE,
                            )
            except IndexError:
                pass

        return G

    # TODO: change 'start' and 'end' to left and right end for regions...

    @staticmethod
    def alignment_hash(a):
        return (
            a.query_region.left_end,
            a.query_region.right_end,
            a.query_region.direction,
            a.type,
        )

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
