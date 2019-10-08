"""Alignments."""
from __future__ import annotations

from collections.abc import Sized
from itertools import count
from typing import List
from typing import Union

from dasi.constants import Constants
from dasi.exceptions import AlignmentException
from dasi.utils import Region


class Alignment(Sized):
    """A pair of Regions that 'aligns' two regions of DNA sequences. All
    regions must always be the same length.

    A subregion of both regions may be taken.
    """

    __slots__ = [
        "query_region",
        "subject_region",
        "type",
        "query_key",
        "subject_key",
        "uid",
    ]

    def __init__(
        self,
        query_region: Region,
        subject_region: Region,
        atype: str,
        query_key: str,
        subject_key: str,
    ):
        """Makes an alignment between two regions of sequences. Validates the
        regions are the same length.

        :param query_region: Query region this alignment aligns to
        :param subject_region: Subject region this alignment aligns to.
        :param atype: Type of alignment
        :param query_key: The record identifier for the query
        :param subject_key: The record identifier for the subject
        """
        self.query_region = query_region
        self.subject_region = subject_region
        self.validate()
        self.type = atype
        self.query_key = query_key
        self.subject_key = subject_key

    def validate(self):
        if not len(self.query_region) == len(self.subject_region):
            raise AlignmentException(
                "Regions must have the same size: {} vs {}. {} vs {}".format(
                    len(self.query_region),
                    len(self.subject_region),
                    self.query_region,
                    self.subject_region,
                )
            )

    def is_perfect_subject(self):
        return len(self.subject_region) == self.subject_region.context_length

    def sub_region(self, qstart: int, qend: int, atype=None) -> Alignment:
        """Returns a copy of the alignment between the inclusive start and end
        relative to the query region.

        :param qstart: start of the query sub region
        :param qend: end of the query sub region
        :param atype: optional type of alignment to return
        :return:
        """
        query_copy = self.query_region.sub(qstart, qend)
        i = self.query_region.i(qstart)
        if i < 0:
            i = self.query_region.i(qstart + self.query_region.context_length)
        if i == len(self.subject_region):
            i = 0
        subject_copy = self.subject_region[i : i + len(query_copy)]

        if atype is None:
            atype = self.type
        self.validate()
        return self.__class__(
            query_region=query_copy,
            subject_region=subject_copy,
            atype=atype,
            query_key=self.query_key,
            subject_key=self.subject_key,
        )

    def copy(self, atype=None) -> Alignment:
        """Do shallow copy of this alignment. Query and subject regions between
        this and the copied alignment will be identical.

        :param atype: new alignment type
        :return:
        """
        if atype is None:
            atype = self.type
        return self.__class__(
            self.query_region,
            self.subject_region,
            atype,
            self.query_key,
            self.subject_region,
        )

    def __len__(self) -> int:
        return len(self.query_region)

    def __str__(self) -> str:
        return "<{} {} {} {}>".format(
            self.__class__.__name__, self.type, self.query_region, self.subject_region
        )

    # def __setstate__(self, state):
    #     self._registry[self.uid] = self

    def __repr__(self) -> str:
        return str(self)


class AlignmentGroupBase:
    """A representative Alignment representing a group of alignments."""

    __slots__ = ["query_region", "alignments", "name", "type", "meta"]

    def __init__(
        self,
        alignments: List[Alignment],
        group_type: str,
        name: str = None,
        query_region: Region = None,
        meta: dict = None,
    ):
        """

        :param alignments:
        :param group_type:
        :param name:
        :param query_region:
        :param meta:
        """
        self.query_region = query_region
        self.alignments = alignments
        self.name = name
        self.type = group_type
        self.meta = meta

    @property
    def query_key(self) -> str:
        """Return the query key associated with the query region."""
        return self.alignments[0].query_key

    @property
    def subject_regions(self) -> List[Region]:
        """Return the list of subject regions in this alignment group."""
        return [a.subject_region for a in self.alignments]

    @property
    def subject_keys(self) -> List[str]:
        """Return the list of subject keys in this alignment group."""
        return [a.subject_key for a in self.alignments]

    def sub_region(self, qstart: int, qend: int, atype: str) -> AlignmentGroupBase:
        """Produce a new alignment group with sub-regions of the query region
        and subject regions at the specified new indicies."""
        alignments_copy = []
        for a in self.alignments:
            alignments_copy.append(a.sub_region(qstart, qend))
        for a in alignments_copy:
            a.type = atype
        return self.__class__(
            alignments=alignments_copy, group_type=atype, name="subregion"
        )

    def __repr__(self) -> str:
        return "<AlignmentGroup {}>".format(self.query_region)


class AlignmentGroup(AlignmentGroupBase):
    """A representative Alignment representing a group of alignments sharing
    the same starting and ending position for a query sequence."""

    __slots__ = ["query_region", "alignments", "name", "type"]

    def __init__(
        self,
        alignments: List[Alignment],
        group_type: str,
        name: str = None,
        meta: dict = None,
    ):
        super().__init__(
            alignments=alignments,
            group_type=group_type,
            name=name,
            query_region=alignments[0].query_region,
            meta=meta,
        )


# class ComplexAlignmentGroup(AlignmentGroupBase):
#     """A representation of an alignment in which the query region is the
#     concatenation of the underlying alignments provided."""
#
#     __slots__ = ["query_region", "alignments", "name", "type", "meta"]
#
#     def __init__(self, alignments: List[Alignment], group_type: str, meta: dict = None):
#         # TODO: adjust alignments
#         query_region = alignments[0].query_region
#         query_region = query_region.new(
#             alignments[0].query_region.a, alignments[-1].query_region.b
#         )
#         super().__init__(
#             alignments=alignments,
#             group_type=group_type,
#             query_region=query_region,
#             meta=meta,
#         )


class PCRProductAlignmentGroup(AlignmentGroupBase):
    """Represents a PCR product alignment from a template alignment and
    forward and reverse alignments. Represents several situations:

    ::

            Situations the PCRProductAlignmentGroup represents:

            Rev primer with overhang
                        <--------
            ------------------
            ---->

            Primers with overhangs
                           <--------
               ------------------
            -------->

            Primers 'within' the template
                    <-----
               ------------------
            -------->

            And so on...

    """

    def __init__(
        self,
        fwd: Union[None, Alignment],
        template: Alignment,
        rev: Union[None, Alignment],
        group_type: str,
        meta: dict = None,
    ):
        """Initialize a new PCRProductAlignmentGroup. Represents a PCR product.
        Query region end points are determined from the first non-None
        alignment and the last non-None alignment. Produces *one new
        alignment*, the template alignment, which is the intersection of the
        provided template alignment and the query_region, which represents the
        exact region for which PCR primers ought to align in a PCR reaction.

        :param fwd: the forward primer alignment. Can be 'inside' the template or have an overhang.
        :param template: template alignment
        :param rev: the reverse primer alignment. Can be 'within' the template or have an overhang.
        :param group_type: group type name
        :param meta: extra meta data
        """
        if fwd is None and rev is None:
            raise AlignmentException("Must provide either a fwd and/or rev alignments")
        alignments = [x for x in [fwd, template, rev] if x is not None]
        a = alignments[0].query_region.a
        b = alignments[-1].query_region.b

        query_region = template.query_region.new(a, b)
        self.raw_template = template
        self._template = None
        self.fwd = fwd
        self.rev = rev

        alignments = [
            x for x in [self.fwd, self.raw_template, self.rev] if x is not None
        ]
        super().__init__(
            alignments=alignments,
            group_type=group_type,
            query_region=query_region,
            meta=meta,
        )

    def get_template(self):
        if self._template is None:
            intersection = self.raw_template.query_region.intersection(
                self.query_region
            )
            self._template = self.raw_template.sub_region(
                intersection.a, intersection.b
            )
        return self._template

    @property
    def subject_keys(self):
        raise AlignmentException(
            "Use subject keys directly, as in `self.fwd.subject_key`"
        )


# TODO: rename this class
class MultiPCRProductAlignmentGroup(AlignmentGroupBase):
    """A PCR Product Alignment with redundant forward primer, reverse primer,
    and template alignments."""

    def __init__(
        self,
        fwds: List[Alignment],
        templates: List[Alignment],
        revs: List[Alignment],
        query_region: Region,
        group_type: str,
    ):
        """Initializes a new MultiPCRProductAlignmentGroup.

        :param fwds:
        :param templates:
        :param revs:
        :param query_region:
        :param group_type:
        """
        self.fwds = fwds
        self.revs = revs
        self.raw_templates = templates
        self._templates = [None] * len(self.raw_templates)
        alignments = [
            x for x in self.fwds + self.raw_templates + self.revs if x is not None
        ]
        super().__init__(
            alignments=alignments, query_region=query_region, group_type=group_type
        )

    def get_template(self, index):
        if self._templates[index] is None:
            intersection = self.raw_templates[index].query_region.intersection(
                self.query_region
            )
            self._templates[index] = self.raw_templates[index].sub_region(
                intersection.a, intersection.b
            )
        return self._templates[index]
