from copy import copy
from typing import Tuple

from numpy import inf


class Constants:
    """DASi constants."""

    FRAGMENT = (
        "PRE-MADE DNA FRAGMENT"
    )  #: an alignment that is generate from an already existing PCR product or fragment
    PCR_PRODUCT = (
        "PCR_PRODUCT"
    )  #: an alignment that is to be generated from a PCR product
    PCR_PRODUCT_WITH_PRIMERS = (
        "PCR_PRODUCT_WITH_PRIMERS"
    )  #: PCR product that can be produces from existing primers
    PCR_PRODUCT_WITH_LEFT_PRIMER = (
        "PCR_PRODUCT_WITH_LEFT_PRIMER"
    )  #: PCR product with existing left primer
    PCR_PRODUCT_WITH_RIGHT_PRIMER = (
        "PCR_PRODUCT_WITH_RIGHT_PRIMER"
    )  #: PCR product with existing right primer
    SHARED_FRAGMENT = (
        "FRAGMENT_SHARED_WITH_OTHER_QUERIES"
    )  #: A fragment alignment that is shared with other queries for potential reuse
    PRIMER_EXTENSION_PRODUCT = (
        "PRIMER_EXTENSION_PRODUCT"
    )  #: a pcr product that uses no template, but extends two primers
    PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER = (
        "PRIMER_EXTENSION_PRODUCT"
    )  #: a pcr product that uses no template, but extends an existing fwd primer.
    PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER = (
        "PRIMER_EXTENSION_PRODUCT"
    )  #: a pcr product that uses no template, but extends an existing rev primer.

    TEMPLATE = "TEMPLATE"  #: A template alignment group. Not an actual molecule.
    GAP = "GAP"  #: region that represents a gap that must be synthesized
    OVERLAP = "OVERLAP"  #: region that represents overlapping molecules
    MISSING = "__MISSING"  #: missing region

    PRIMER = "PRIMER"  #: a primer binding alignment

    PRIMER_MIN_BIND = 14  #: minimum primer binding for searching for primer alignments
    MIN_OVERLAP = 15  #: minimum overlap for searching for overlapping alignments
    MAX_HOMOLOGY = 100  #: maximum overlap for searching for overlapping alignments
    INF = 10.0 ** 6  #: almost infinity


class MoleculeType:

    types = {}

    def __init__(
        self,
        name: str,
        design: Tuple[bool, bool],
        use_direct: bool,
        cost: float,
        efficiency=1.0,
        synthesize=False,
    ):
        self.name = name
        self.cost = cost
        self.use_direct = use_direct
        self.synthesize = synthesize
        self.types[name] = self
        self.efficiency = efficiency
        self.design = design
        self.int_or_ext = None

    def __repr__(self):
        return "<{} name='{}'>".format(self.__class__.__name__, self.name)


class InternalType(MoleculeType):
    def __init__(
        self,
        name,
        design,
        use_direct: bool,
        cost: float,
        efficiency=1.0,
        synthesize: bool = False,
    ):
        super().__init__(
            name,
            use_direct=use_direct,
            design=design,
            cost=cost,
            efficiency=efficiency,
            synthesize=synthesize,
        )
        self.design = design
        self.int_or_ext = "internal"


class ExternalType(MoleculeType):
    def __init__(
        self,
        name,
        use_direct: bool,
        cost: float,
        efficiency=1.0,
        synthesize: bool = False,
    ):
        super().__init__(
            name,
            use_direct=use_direct,
            design=None,
            cost=cost,
            efficiency=efficiency,
            synthesize=synthesize,
        )
        self.int_or_ext = "external"

    def __call__(self, design):
        copied = copy(self)
        copied.design = design
        return copied


InternalType(Constants.FRAGMENT, (False, False), True, 0.0, 0.98)
InternalType(Constants.PCR_PRODUCT, (True, True), False, 10.0, 0.95)
InternalType(Constants.PCR_PRODUCT_WITH_PRIMERS, (False, False), False, 10.0, 0.95)
InternalType(Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, (True, False), False, 10.0, 0.95)
InternalType(Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, (False, True), False, 10.0, 0.95)
InternalType(Constants.PRIMER_EXTENSION_PRODUCT, (False, False), False, 9.0, 0.95)
InternalType(
    Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER, (False, True), False, 9.0, 0.95
)
InternalType(
    Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
    (True, False),
    False,
    9.0,
    0.95,
)
# TODO: implement primer dimer product with left and right primers

ExternalType(
    name=Constants.OVERLAP, use_direct=False, cost=0.0, efficiency=1.0, synthesize=False
)
ExternalType(
    name=Constants.GAP, use_direct=False, cost=0.0, efficiency=1.0, synthesize=True
)

MoleculeType(
    name=Constants.MISSING,
    design=None,
    use_direct=False,
    cost=inf,
    efficiency=0.0,
    synthesize=False,
)


MoleculeType(
    name=Constants.PRIMER,
    design=None,
    use_direct=False,
    cost=inf,
    efficiency=0.0,
    synthesize=False,
)
