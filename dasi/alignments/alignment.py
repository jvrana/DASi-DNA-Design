"""
Alignments
"""

from dasi.utils import Region
from dasi.exceptions import AlignmentException
from typing import List
from collections.abc import Sized

ALIGNMENT_SLOTS = [
    "query_region",
    "subject_region",
    "type",
    "query_key",
    "subject_key",
    "grouping_tags",
]


class Alignment(Sized):
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
        atype: str,
        query_key: str,
        subject_key: str,
    ):
        """
        Makes an alignment between two regions of sequences. Validates the regions are the same length.

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
        self.grouping_tags = {}

    def validate(self):
        if not len(self.query_region) == len(self.subject_region):
            raise AlignmentException(
                "Regions must have the same size: {} vs {}. {} vs {}".format(
                    len(self.query_region), len(self.subject_region), self.query_region, self.subject_region
                )
            )

    def is_perfect_subject(self):
        return len(self.subject_region) == self.subject_region.bontext_length

    def sub_region(self, qstart: int, qend: int, atype=None):
        """
        Returns a copy of the alignment between the inclusive start and end relative to the
        query region.

        :param qstart: start of the query sub region
        :param qend: end of the query sub region
        :param atype: optional type of alignment to return
        :return:
        """
        query_copy = self.query_region.sub(qstart, qend)
        i = self.query_region.i(qstart)
        if i < 0:
            i = self.query_region.i(qstart + self.query_region.context_length)
        try:
            subject_copy = self.subject_region[i:i+len(query_copy)]
        except Exception as e:
            raise e

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

    def copy(self, atype=None):
        if atype is None:
            atype = self.type
        return self.__class__(
            self.query_region,
            self.subject_region,
            atype,
            self.query_key,
            self.subject_region,
        )

    def __len__(self):
        return len(self.query_region)

    def __str__(self):
        return "<{} {} {} {}>".format(
            self.__class__.__name__, self.type, self.query_region, self.subject_region
        )

    def __repr__(self):
        return str(self)


class AlignmentGroupBase(object):
    """
    A representative Alignment representing a group of alignments.
    """

    __slots__ = ["query_region", "alignments", "name", "type"]

    def __init__(
        self,
        alignments: List[Alignment],
        group_type: str,
        name=None,
        query_region=None
    ):
        self.query_region = query_region
        self.alignments = alignments
        self.name = name
        self.type = group_type

    @property
    def query_key(self):
        return self.alignments[0].query_key

    @property
    def subject_regions(self):
        return [a.subject_region for a in self.alignments]

    @property
    def subject_keys(self):
        return [a.subject_key for a in self.alignments]

    def sub_region(self, qstart: int, qend: int, atype: str):
        alignments_copy = []
        for a in self.alignments:
            alignments_copy.append(a.sub_region(qstart, qend))
        for a in alignments_copy:
            a.type = atype
        return self.__class__(
            alignments=alignments_copy, group_type=atype, name="subregion"
        )

    def __repr__(self):
        return "<AlignmentGroup {}>".format(self.query_region)


class AlignmentGroup(AlignmentGroupBase):
    """
    A representative Alignment representing a group of alignments sharing the
    same starting and ending position for a query sequence.
    """

    __slots__ = ["query_region", "alignments", "name", "type"]

    def __init__(self, alignments: List[Alignment], group_type: str, name=None):
        super().__init__(alignments=alignments, group_type=group_type, name=name, query_region=alignments[0].query_region)


class ComplexAlignmentGroup(AlignmentGroupBase):
    """
    A representation of an alignment in which the query region is the concatenation of the
    underlying alignments provided.
    """

    __slots__ = ["query_region", "alignments", "name", "type"]

    def __init__(self, alignments: List[Alignment], group_type: str):
        query_region = alignments[0].query_region
        query_region = query_region.new(
            alignments[0].query_region.a, alignments[-1].query_region.b
        )
        super().__init__(alignments=alignments, group_type=group_type, query_region=query_region)
