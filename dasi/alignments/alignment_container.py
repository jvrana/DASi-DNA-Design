from .alignment import Alignment, AlignmentGroup, ComplexAlignmentGroup
from dasi.log import logger
from dasi.utils import Region, bisect_slice_between, sort_with_keys
from dasi.constants import Constants
from dasi.exceptions import AlignmentContainerException
from more_itertools import partition, unique_everseen, flatten
from typing import Dict, List
from Bio.SeqRecord import SeqRecord
from bisect import bisect_left
from copy import deepcopy
from collections.abc import Sized
from uuid import uuid4

def blast_to_region(query_or_subject, seqdb):
    """
    Converts a blast data result to a Region. Blast
    results are indicated by two positions with index starting
    at 1 and positions being inclusive. This returns a Region
    with index starting at 0 and the end point position being
    exclusive.

    :param query_or_subject:
    :type query_or_subject:
    :param seqdb:
    :type seqdb:
    :return:
    :rtype:
    """
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

    def _create_pcr_product_alignment(self, template_group: AlignmentGroup, fwd: Alignment, rev: Alignment, alignment_type: str):
        groups = []
        for a in template_group.alignments:
            alignments = [x for x in [fwd, a, rev] if x is not None]
            self._new_grouping_tag(alignments, alignment_type)
            groups.append(ComplexAlignmentGroup(alignments, alignment_type))
        return groups

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
        fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.b)
        rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.a)
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
                    try:
                        _rev_span = g.query_region.sub(f.query_region.a, g.query_region.b)
                    except IndexError:
                        _rev_span = deepcopy(g.query_region)
                    for a, b in _rev_span.ranges():
                        _rev_bind += bisect_slice_between(
                            rev_bind, rkeys, a, b
                        )
                for r in _rev_bind:
                    pairs += self._create_pcr_product_alignment(g, f, r, Constants.PCR_PRODUCT_WITH_PRIMERS)

            # left primer
            for f in fwd_bind:
                pairs += self._create_pcr_product_alignment(g, f, None, Constants.PCR_PRODUCT_WITH_LEFT_PRIMER)

            # right primer
            for r in rev_bind:
                pairs += self._create_pcr_product_alignment(g, None, r, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER)
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

            overlapping = self.filter_alignments_by_span(
                group_sort,
                group_a.query_region,
                key=lambda p: p.query_region.a)

            #
            for group_b in overlapping:
                if group_b is not group_a:
                    left = group_a.sub_region(group_a.query_region.a, group_b.query_region.a, type)
                    overlap = group_a.sub_region(group_b.query_region.a, group_a.query_region.b, type)
                    right = group_b.sub_region(group_a.query_region.b, group_b.query_region.b, type)

                    if len(left.query_region) > MIN_OVERLAP:
                        alignments += left.alignments

                    if len(overlap.query_region) > MIN_OVERLAP:
                        alignments += overlap.alignments

                    if len(right.query_region) > MIN_OVERLAP:
                        alignments += right.alignments
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
            "Number of total groups: {}".format(len(self.groups()))
        )

    @classmethod
    def _new_grouping_tag(cls, alignments, type: str, key=None):
        """
        Make a new ordered grouping by type and a uuid.

        :param alignments:
        :type alignments:
        :param type:
        :type type:
        :param key:
        :type key:
        :return:
        :rtype:
        """
        if key is None:
            key = str(uuid4())
        group_key = (key, type)
        for i, a in enumerate(alignments):
            if key in a.grouping_tags:
                raise AlignmentContainerException("Key '{}' already exists in grouping tag".format(key))
            a.grouping_tags[group_key] = (i)


    # TODO: change 'start' and 'end' to left and right end for regions...
    # TODO: type increases computation time exponentially, we only need 'b' from group1 and 'a' from group2
    @staticmethod
    def _alignment_hash(a):
        """A hashable representation of an alignment for grouping."""
        return (a.query_region.a, a.query_region.b, a.query_region.direction, a.type)

    @classmethod
    def complex_alignment_groups(cls, alignments: List[Alignment]) -> List[AlignmentGroup]:
        key_to_alignments = {}
        for a in alignments:
            for (uuid, type), i in a.grouping_tags.items():
                key_to_alignments.setdefault((uuid, type), list()).append((i, a))
        complex_groups = []
        for (uuid, type), alist in key_to_alignments:
            sorted_alist = [x[-1] for x in sorted(alist)]
            complex_groups.append(ComplexAlignmentGroup(sorted_alist, type))
        return complex_groups

    @classmethod
    def redundent_alignment_groups(cls, alignments: List[Alignment]) -> List[AlignmentGroup]:
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
            alignment_groups.append(AlignmentGroup(group, group[0].type))
        return alignment_groups

    def groups(self) -> List[AlignmentGroup]:
        allgroups = []
        allgroups += self.redundent_alignment_groups(self.alignments)
        allgroups += self.complex_alignment_groups(self.alignments)
        return allgroups

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

    # TODO: change from property
    @property
    def groups_by_type(self) -> Dict[str, AlignmentGroup]:
        """
        Return alignment groups according to their alignment 'type'

        :return: dict
        """
        d = {}
        for g in self.groups():
            d.setdefault(g.type, list()).append(g)
        return d

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


