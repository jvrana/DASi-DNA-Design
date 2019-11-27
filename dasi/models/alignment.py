"""Alignments."""
from collections.abc import Sized
from typing import Dict
from typing import List
from typing import Union

from .molecule import MoleculeType
from dasi.exceptions import AlignmentException
from dasi.utils import Region


class RepresentsMolecule:
    """Mixin for molecular sized alignments or alignment groups."""

    def __init__(self, query_region: Region, atype: str):
        self.query_region = query_region
        self.type = atype
        if atype not in MoleculeType.types:
            raise ValueError("atype '{}' not in MoleculeTypes".format(atype))

    def size_ok(self):
        size = len(self.query_region)
        mol_type = MoleculeType.types[self.type]
        if mol_type.min_size is not None and size < mol_type.min_size:
            return False
        if mol_type.max_size is not None and size > mol_type.max_size:
            return False
        return True


class Alignment(RepresentsMolecule, Sized):
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
        meta: dict = None,
    ):
        """Makes an alignment between two regions of sequences. Validates the
        regions are the same length.

        :param query_region: Query region this alignment aligns to
        :param subject_region: Subject region this alignment aligns to.
        :param atype: Type of alignment
        :param query_key: The record identifier for the query
        :param subject_key: The record identifier for the subject
        """
        super().__init__(query_region, atype)
        assert query_region.direction == 1
        self.subject_region = subject_region
        self.validate()
        self.query_key = query_key
        self.subject_key = subject_key
        if meta is None:
            meta = {}
        self.meta = meta

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

    def sub_region(self, qstart: int, qend: int, atype=None) -> "Alignment":
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

        if self.subject_region.direction == -1:
            b = len(self.subject_region) - i
            a = b - len(query_copy)
            if a == b == len(self.subject_region):
                a = b = 0
            subject_copy = self.subject_region[a:b]
        else:
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

    def copy(self, atype=None) -> "Alignment":
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
            self.subject_key,
        )

    def __len__(self) -> int:
        return len(self.query_region)

    def __str__(self) -> str:
        return "<{} {} {} {}>".format(
            self.__class__.__name__, self.type, self.query_region, self.subject_region
        )

    @staticmethod
    def _rhash(region: Region):
        return (region.a, region.b, region.c, region.context_length, region.cyclic)

    def eq_hash(self):
        return (
            (self.query_key, self.subject_key, self.type)
            + self._rhash(self.query_region)
            + self._rhash(self.subject_region)
        )

    def __eq__(self, other: "Alignment"):
        return self.eq_hash() == other.eq_hash()

    def __repr__(self) -> str:
        return str(self)


class AlignmentGroupBase(RepresentsMolecule):
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
        super().__init__(query_region, group_type)
        self.alignments = alignments
        self.name = name
        if meta is None:
            meta = {}
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

    def sub_region(self, qstart: int, qend: int, atype: str) -> "AlignmentGroupBase":
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
        query_region: Region,
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

        query_region = query_region.new(a, b)

        self.raw_template = template
        self._template = None
        self.fwd = fwd
        self.rev = rev

        super().__init__(
            alignments=alignments,
            group_type=group_type,
            query_region=query_region,
            meta=meta,
        )

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
        groupings: List[Dict[str, Alignment]],
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
        self.groupings = groupings
        self._templates = [None] * len(self.groupings)
        fwds = [d["fwd"] for d in self.groupings]
        revs = [d["rev"] for d in self.groupings]
        templates = [d["template"] for d in self.groupings]
        alignments = [a for a in fwds + revs + templates if a is not None]
        super().__init__(
            alignments=alignments, query_region=query_region, group_type=group_type
        )

    def get_template(self, index):
        if self._templates[index] is None:
            template = self.groupings[index]["template"]
            intersection = template.query_region.intersection(self.query_region)
            self._templates[index] = template.sub_region(intersection.a, intersection.b)
        return self._templates[index]
