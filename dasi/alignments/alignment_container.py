"""Alignment container."""
from bisect import bisect_left
from collections.abc import Sized
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

from Bio.SeqRecord import SeqRecord
from frozendict import frozendict
from more_itertools import flatten
from more_itertools import partition
from more_itertools import unique_everseen

from .alignment import Alignment
from .alignment import AlignmentGroup
from .alignment import PCRProductAlignmentGroup
from dasi.constants import Constants
from dasi.constants import MoleculeType
from dasi.exceptions import AlignmentContainerException
from dasi.log import logger
from dasi.utils import bisect_slice_between
from dasi.utils import Region
from dasi.utils import sort_with_keys


def blast_to_region(query_or_subject, seqdb):
    """Converts a blast data result to a Region. Blast results are indicated by
    two positions with index starting at 1 and positions being inclusive. This
    returns a Region with index starting at 0 and the end point position being
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

    s, e = data["start"], data["raw_end"]
    length = len(record)

    if data["strand"] == -1:
        s, e = e, s

    region = Region(
        s - 1,
        e,
        length=length,
        cyclic=data["circular"],
        direction=data["strand"],
        index=0,
        name="{}: {}".format(record.id, record.name),
        ignore_wrap=False,
    )
    return region


class AlignmentContainer(Sized):
    """Container for a set of query-to-subject alignments for a single query.

    Instance Attributes/Properties:
        alignments  list of alignments for the container
                    List[Alignment]
        seqdb       key to SeqRecord dictionary
                    Dict[str, SeqRecord]
        logger      the instance's logger
                    Loggable

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
        Constants.SHARED_FRAGMENT,
        Constants.GAP,
        Constants.OVERLAP,
        Constants.MISSING,
    )  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord], alignments=None):
        self._alignments = []
        self._frozen = False
        if alignments is not None:
            self.alignments = alignments
        self.seqdb = seqdb
        self.logger = logger(self)

    @property
    def alignments(self):
        return self._alignments

    @alignments.setter
    def alignments(self, v):
        if self._frozen:
            raise AlignmentContainerException(
                "Cannot set alignments. Container is frozen."
            )
        self._alignments = v
        self._check_single_query_key(self._alignments)

    @staticmethod
    def _check_single_query_key(alignments):
        keys = {a.query_key for a in alignments}
        if len(keys) > 1:
            raise AlignmentContainerException(
                "AlignmentContainer cannot contain more than one query. Contains the "
                "following query keys: {}".format(keys)
            )

    @classmethod
    def filter_alignments_by_span(
        cls, alignments, region, key=None, end_inclusive=True
    ):
        fwd, fwd_keys = sort_with_keys(alignments, key=key)
        found = []
        for a, b in region.ranges(ignore_wraps=True):
            if not end_inclusive:
                b = b - 1
            found += bisect_slice_between(fwd, fwd_keys, a, b)
        return found

    def find_groups_by_pos(
        self, a: int, b: int, group_type: str, groups=None
    ) -> List[Union[AlignmentGroup, PCRProductAlignmentGroup]]:
        """Return a list of groups that have the same positions (a, b) and same
        group_type.

        :param a: starting position
        :param b: ending position
        :param group_type: group_type
        :param groups: optional list of groups to search
        :return: list of groups
        """

        if group_type not in list(self.valid_types) + ["ANY"]:
            raise AlignmentContainerException(
                "Type '{}' not found in valid types: {}".format(
                    group_type, self.valid_types
                )
            )
        if groups is None:
            groups = self.groups()
        found = []
        for g in groups:
            if (
                g.query_region.a == a
                and g.query_region.b == b
                and (group_type == "ANY" or g.type == group_type)
            ):
                found.append(g)
        return found

    def _create_pcr_product_alignment(
        self,
        template_group: AlignmentGroup,
        fwd: Union[Alignment, None],
        rev: Union[Alignment, None],
        alignment_type: str,
    ):
        """
        Create a new alignment group for a PCR product.

        ::

            Situations the PCRProductAlignmentGroup represents:

                        <--------
            ------------------
            ---->

                           <--------
               ------------------
            -------->


                    <-----
               ------------------
            -------->


                          <------
               ------------------
            -------->

        :param template_group:
        :param fwd:
        :param rev:
        :param alignment_type:
        :return:
        """
        groups = []
        for a in template_group.alignments:
            # TODO: here, you should take a subregion and set it as a new type
            #       that is not evaluated in other contexts
            product_group = PCRProductAlignmentGroup(
                fwd=fwd, template=a, rev=rev, group_type=alignment_type
            )
            alignments = (
                product_group.fwd,
                product_group.raw_template,
                product_group.rev,
            )
            self._new_grouping_tag(alignments, alignment_type)
            self.alignments += product_group.alignments
        return groups

    def expand_primer_pairs(
        self, alignment_groups: List[AlignmentGroup]
    ) -> List[Alignment]:
        """Creates new alignments for all possible primer pairs. Searches for
        fwd and rev primer pairs that exist within other alignments and
        produces all combinations of alignments that can form from these primer
        pairs.

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
            fwd_bind_region = g.query_region[Constants.PRIMER_MIN_BIND :]
            rev_bind_region = g.query_region[: -Constants.PRIMER_MIN_BIND]
            fwd_bind = self.filter_alignments_by_span(
                fwd, fwd_bind_region, key=lambda p: p.query_region.b
            )
            rev_bind = self.filter_alignments_by_span(
                rev, rev_bind_region, key=lambda p: p.query_region.a
            )
            rev_bind, rkeys = sort_with_keys(rev_bind, key=lambda p: p.query_region.a)

            # both primers
            for f in fwd_bind:
                _rev_bind = []
                if len(query_ranges) == 1:
                    i = bisect_left(rkeys, f.query_region.a)
                    _rev_bind = rev_bind[i:]
                else:
                    try:
                        _rev_span = g.query_region.sub(
                            f.query_region.a, g.query_region.b
                        )
                    except IndexError:
                        _rev_span = g.query_region

                    for a, b in _rev_span.ranges():
                        _rev_bind += bisect_slice_between(rev_bind, rkeys, a, b)
                for r in _rev_bind:
                    pairs += self._create_pcr_product_alignment(
                        g, f, r, Constants.PCR_PRODUCT_WITH_PRIMERS
                    )

            # left primer
            for f in fwd_bind:
                pairs += self._create_pcr_product_alignment(
                    g, f, None, Constants.PCR_PRODUCT_WITH_LEFT_PRIMER
                )

            # right primer
            for r in rev_bind:
                pairs += self._create_pcr_product_alignment(
                    g, None, r, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER
                )
        return pairs

    def expand_overlaps(
        self, alignment_groups: List[AlignmentGroup], atype=Constants.PCR_PRODUCT
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

        :param alignment_groups: list of alignment groups to expand
        :param atype: the alignment type label for expanded alignments
        :return: list
        """

        min_overlap = Constants.MIN_OVERLAP
        group_sort, group_keys = sort_with_keys(
            alignment_groups, key=lambda x: x.query_region.a
        )
        alignments = []
        for group_a in logger.tqdm(group_sort, "INFO", desc="expanding pcr products"):
            overlapping = self.filter_alignments_by_span(
                group_sort,
                group_a.query_region,
                key=lambda p: p.query_region.a,
                end_inclusive=False,
            )

            for group_b in overlapping:
                if group_b is not group_a:
                    left = group_a.sub_region(
                        group_a.query_region.a, group_b.query_region.a, atype
                    )
                    overlap = group_a.sub_region(
                        group_b.query_region.a, group_a.query_region.b, atype
                    )

                    if len(left.query_region) > min_overlap:
                        alignments += left.alignments

                    if len(overlap.query_region) > min_overlap:
                        alignments += overlap.alignments
        return alignments

    # TODO: break apart long alignments
    def expand(self, expand_overlaps=True, expand_primers=True):
        """Expand the number of alignments in this container using overlaps or
        primers.

        :param expand_overlaps:
        :param expand_primers:
        :return:
        """

        self.logger.info("=== Expanding alignments ===")
        # We annotate any original PCR_PRODUCT with FRAGMENT if they are
        # 'perfect_subjects'
        # This means they already exist as pre-made fragments
        self.logger.info("Number of alignments: {}".format(len(self.alignments)))

        templates = self.get_groups_by_types(
            [Constants.PCR_PRODUCT, Constants.FRAGMENT]
        )

        if expand_primers:
            pairs = self.expand_primer_pairs(templates)
            self.logger.info("Number of pairs: {}".format(len(pairs)))

        if expand_overlaps:
            expanded = self.expand_overlaps(templates)
            self.alignments += expanded
            self.logger.info("Number of new alignments: {}".format(len(expanded)))
        self.logger.info("Number of total alignments: {}".format(len(self.alignments)))
        self.logger.info("Number of total groups: {}".format(len(self.groups())))

    @classmethod
    def _new_grouping_tag(
        cls, alignments: List[Alignment], atype: str, key: Any = None, meta: dict = None
    ):
        """Make a new ordered grouping by type and a uuid.

        Grouping tags are maintained as an attribute on the Alignment class.
        The grouping tags are of the form:
            {(<uuid>, <group_type>): (<index>, <metadata>)}
        which registers an Alignment to a specific group in the container.
        `None` values in the alignments are skipped.

        :param alignments: list of alignments
        :param atype:
        :type atype:
        :param key: optional key
        :return:
        :rtype:
        """
        if key is None:
            key = str(uuid4())
        group_key = (key, atype)
        for i, a in enumerate(alignments):
            if a is not None:
                if key in a.grouping_tags:
                    raise AlignmentContainerException(
                        "Key '{}' already exists in grouping tag".format(key)
                    )
                a.grouping_tags[group_key] = i, meta

    @staticmethod
    def _alignment_hash(a):
        """A hashable representation of an alignment for grouping."""
        return a.query_region.a, a.query_region.b, a.query_region.direction, a.type

    @classmethod
    def pcr_alignment_groups(
        cls, alignments: List[Alignment]
    ) -> List[Union[AlignmentGroup, PCRProductAlignmentGroup]]:
        key_to_alignments = {}
        for a in alignments:
            if not isinstance(a, Alignment):
                raise Exception
        for a in alignments:
            for (uuid, atype), (i, meta) in a.grouping_tags.items():
                key_to_alignments.setdefault((uuid, atype), list()).append((i, a))
        complex_groups = []
        for (uuid, atype), alist in key_to_alignments.items():
            alist_dict = dict(x for x in alist)
            complex_groups.append(
                PCRProductAlignmentGroup(
                    fwd=alist_dict.get(0, None),
                    template=alist_dict[1],
                    rev=alist_dict.get(2, None),
                    group_type=atype,
                )
            )
        return complex_groups

    @classmethod
    def redundent_alignment_groups(
        cls, alignments: List[Alignment]
    ) -> List[AlignmentGroup]:
        """Return AlignmentGroups that have been grouped by alignment_hash.

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
        allgroups += self.pcr_alignment_groups(self.alignments)
        return allgroups

    @property
    def types(self) -> Tuple[Any, ...]:
        """Return all valid types.

        :return:
        """
        return tuple(self.valid_types)

    def get_groups_by_types(
        self, types: List[str]
    ) -> Union[AlignmentGroup, List[AlignmentGroup]]:
        """Return AlignmentGroups by fragment type.

        :param types: list of types
        :return:
        """
        groups = self.groups_by_type()
        if isinstance(types, str):
            return groups[types]
        else:
            return list(unique_everseen(flatten([groups[t] for t in types])))

    def get_alignments_by_types(self, types: List[str]) -> List[Alignment]:
        groups = self.get_groups_by_types(types)
        return list(flatten([g.alignments for g in groups]))

    def groups_by_type(
        self
    ) -> Dict[str, List[Union[AlignmentGroup, PCRProductAlignmentGroup]]]:
        """Return alignment groups according to their alignment 'type'.

        :return: dict
        """
        d = {}
        for t in self.valid_types:
            d[t] = []
        for g in self.groups():
            d[g.type].append(g)
        return d

    def freeze(self):
        """Freeze the container, disallowing further modifications to
        alignments."""
        self._alignments = tuple(self._alignments)
        self._frozen = True

    def unfreeze(self):
        """Unfreeze the container, allowing modifications to alignments."""
        self._alignments = list(self._alignments)
        self._frozen = False

    def __len__(self):
        return len(self.alignments)


class AlignmentContainerFactory:
    """Class that maintains a shared list of alignments and shared sequence
    database.

    AlignmentContainers can be retrieved in a dict grouped by their
    query via `.containers()`
    """

    valid_types = (
        Constants.PCR_PRODUCT,
        Constants.PRIMER,
        Constants.FRAGMENT,
        Constants.PCR_PRODUCT_WITH_PRIMERS,
        Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
        Constants.SHARED_FRAGMENT,
    )  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord]):
        """Construct a new AlignmentContainer.

        :param seqdb: a sequence record database
        """
        self._alignments = (
            {}
        )  # dictionary of query_key to alignment; Dict[str, List[Alignment]]
        self._containers = None
        self.logger = logger(self)
        self.seqdb = seqdb

    @property
    def alignments(self) -> frozendict:
        """Return dict of alignments keyed by query_key.

        :return:
        """
        return frozendict(self._alignments)

    def set_alignments(self, alignments: Dict[str, List[Alignment]]) -> None:
        """Set the alignments.

        :param alignments: new iterable of alignments
        :return:
        """
        self._alignments = alignments
        self._containers = None

    def load_blast_json(self, data: List[Dict], atype: str):
        """Create alignments from a formatted BLAST JSON result.

        :param data: formatted BLAST JSON result
        :param atype: the type of alignment to initialize
        :return: None
        """
        self.logger.info(
            'Loading blast json ({} entries) to fragment type "{}"'.format(
                len(data), atype
            )
        )
        assert atype in self.valid_types
        for d in data:
            query_region = blast_to_region(d["query"], self.seqdb)
            subject_region = blast_to_region(d["subject"], self.seqdb)
            query_key = d["query"]["origin_key"]
            subject_key = d["subject"]["origin_key"]

            alignment = Alignment(
                query_region,
                subject_region,
                atype=atype,
                query_key=query_key,
                subject_key=subject_key,
            )
            self._alignments.setdefault(query_key, list()).append(alignment)

    def containers(self) -> Dict[str, AlignmentContainer]:
        """Return dictionary of AlignmentContainers keyed by query_keys.

        :return:
        """
        if self._containers is None:
            container_dict = {}
            for key, alignments in self.alignments.items():
                container_dict[key] = AlignmentContainer(
                    self.seqdb, alignments=alignments
                )
            self._containers = container_dict
        return frozendict(self._containers)
