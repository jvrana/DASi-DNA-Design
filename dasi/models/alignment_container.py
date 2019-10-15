"""Alignment container."""
from __future__ import annotations

from bisect import bisect_left
from collections.abc import Sized
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

from Bio.SeqRecord import SeqRecord
from frozendict import frozendict
from more_itertools import flatten
from more_itertools import partition
from more_itertools import unique_everseen

from .alignment import Alignment
from .alignment import AlignmentGroup
from .alignment import MultiPCRProductAlignmentGroup
from .alignment import PCRProductAlignmentGroup
from dasi.constants import Constants
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
        Constants.PRIMER_EXTENSION_PRODUCT,
        Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
        Constants.SHARED_FRAGMENT,
        Constants.GAP,
        Constants.OVERLAP,
        Constants.MISSING,
    )  # valid fragment types

    def __init__(self, seqdb: Dict[str, SeqRecord], alignments=None):
        self._alignments = []
        self._frozen = False
        self._frozen_groups = None
        if alignments is not None:
            self.alignments = alignments
        self.seqdb = seqdb
        self.logger = logger(self)
        self.multi_grouping_tags = {}

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
        lim_size: bool,
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
            product_group = PCRProductAlignmentGroup(
                fwd=fwd,
                template=a,
                rev=rev,
                query_region=a.query_region,
                group_type=alignment_type,
            )

            if lim_size and not product_group.size_ok():
                continue
            self._new_multi_pcr_grouping_tag(product_group)
            groups.append(product_group)
        return groups

    def _create_primer_extension_alignment(
        self, fwd: Alignment, rev: Alignment, alignment_type: str, lim_size: bool
    ):
        if fwd is None:
            query_region = rev.query_region
        else:
            query_region = fwd.query_region
        product_group = PCRProductAlignmentGroup(
            fwd=fwd,
            template=None,
            rev=rev,
            query_region=query_region,
            group_type=alignment_type,
        )
        if lim_size and not product_group.size_ok():
            return []
        self._new_multi_pcr_grouping_tag(product_group)
        return [product_group]

    def expand_primer_extension_products(self, only_one_required=False, lim_size=True):
        primers = self.get_alignments_by_types(Constants.PRIMER)

        rev, fwd = partition(lambda p: p.subject_region.direction == 1, primers)
        fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.b)
        rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.a)
        pairs = []
        for f in fwd:
            rev_bind_region = f.query_region[: -Constants.PRIMER_MIN_BIND]
            rev_bind = self.filter_alignments_by_span(
                rev, rev_bind_region, key=lambda p: p.query_region.a
            )
            rev_bind, rkeys = sort_with_keys(rev_bind, key=lambda p: p.query_region.a)

            for r in rev_bind:
                if r.query_region.b in f.query_region:
                    if r.query_region.b == f.query_region.b:
                        pass
                    else:
                        continue
                pairs += self._create_primer_extension_alignment(
                    f, r, Constants.PRIMER_EXTENSION_PRODUCT, lim_size=lim_size
                )

        if only_one_required:
            for f in fwd:
                # existing fwd primer
                pairs += self._create_primer_extension_alignment(
                    f,
                    None,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER,
                    lim_size=lim_size,
                )

            for r in rev:
                # existing fwd primer
                pairs += self._create_primer_extension_alignment(
                    None,
                    r,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
                    lim_size=lim_size,
                )
        return pairs

    def expand_primer_pairs(
        self, alignment_groups: List[AlignmentGroup], lim_size: bool = True
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
                    if f.query_region.a in r.query_region:
                        if f.query_region.a == r.query_region.a:
                            pass
                        else:
                            continue
                    if r.query_region.b in f.query_region:
                        if r.query_region.b == f.query_region.b:
                            pass
                        else:
                            continue
                    pairs += self._create_pcr_product_alignment(
                        g, f, r, Constants.PCR_PRODUCT_WITH_PRIMERS, lim_size=lim_size
                    )
            # left primer
            for f in fwd_bind:
                pairs += self._create_pcr_product_alignment(
                    g,
                    f,
                    None,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                    lim_size=lim_size,
                )

            # right primer
            for r in rev_bind:
                pairs += self._create_pcr_product_alignment(
                    g,
                    None,
                    r,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                    lim_size=lim_size,
                )
        return pairs

    def expand_overlaps(
        self,
        alignment_groups: List[AlignmentGroup],
        atype=Constants.PCR_PRODUCT,
        lim_size: bool = True,
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
        if lim_size:
            alignments = [a for a in alignments if a.size_ok()]
        return alignments

    # TODO: break apart long alignments
    def expand(
        self,
        expand_overlaps=True,
        expand_primers=True,
        expand_primer_dimers=False,
        lim_size: bool = True,
    ):
        """Expand the number of alignments in this container using overlaps or
        primers.

        :param expand_overlaps:
        :param expand_primers:
        :return:
        """

        templates = self.get_groups_by_types(
            [Constants.PCR_PRODUCT, Constants.FRAGMENT]
        )
        if lim_size:
            templates = [t for t in templates if t.size_ok()]

        if expand_primers:
            self.expand_primer_pairs(templates, lim_size=lim_size)

        if expand_primer_dimers:
            self.expand_primer_extension_products(lim_size=lim_size)

        # TODO: why not expand overlaps using the primer pairs???
        if expand_overlaps:
            expanded = self.expand_overlaps(
                templates, atype=Constants.PCR_PRODUCT, lim_size=lim_size
            )
            self.alignments += expanded

        # TODO: make a 'set' of all alignments

    def _new_multi_pcr_grouping_tag(self, group: PCRProductAlignmentGroup):
        group_key = (group.query_region.a, group.query_region.b, group.type)
        self.multi_grouping_tags.setdefault(group_key, list())
        self.multi_grouping_tags[group_key].append(
            {
                "fwd": group.fwd,
                "template": group.raw_template,
                "rev": group.rev,
                "query_region": group.query_region,
            }
        )

    @staticmethod
    def _alignment_hash(a):
        """A hashable representation of an alignment for grouping."""
        return a.query_region.a, a.query_region.b, a.query_region.direction, a.type

    def pcr_alignment_groups(self):
        groups = []
        for (a, b, group_type), adict_list in self.multi_grouping_tags.items():
            g = [
                {"fwd": d["fwd"], "rev": d["rev"], "template": d["template"]}
                for d in adict_list
            ]

            regions = [d["query_region"] for d in adict_list]
            assert len({(q.a, q.b) for q in regions}) == 1

            groups.append(
                MultiPCRProductAlignmentGroup(
                    g, query_region=adict_list[0]["query_region"], group_type=group_type
                )
            )
        return groups

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

    def groups(self) -> List[Union[AlignmentGroup, MultiPCRProductAlignmentGroup]]:
        if self._frozen:
            return self._frozen_groups
        else:
            allgroups = []
            allgroups += self.redundent_alignment_groups(self.alignments)
            allgroups += self.pcr_alignment_groups()
            self._frozen_groups = allgroups
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
        Constants.PRIMER_EXTENSION_PRODUCT,
        Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER,
        Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
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

    def set_alignments(self, alignments):
        self._alignments = alignments
        self._containers = None
