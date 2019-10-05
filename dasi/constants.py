from typing import Tuple
from typing import Union


class Constants:
    """DASi constants."""

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
    SHARED_FRAGMENT = (
        "FRAGMENT_SHARED_WITH_OTHER_QUERIES"
    )  # A fragment alignment that is shared with other queries for potential reuse

    GAP = "GAP"
    OVERLAP = "OVERLAP"
    MISSING = "__MISSING"

    PRIMER = "PRIMER"  # a primer binding alignment

    PRIMER_MIN_BIND = 15
    MIN_OVERLAP = 20
    MAX_HOMOLOGY = 100
    INF = 10.0 ** 6


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
        self.design = design
        self.cost = cost
        self.use_direct = use_direct
        self.synthesize = synthesize
        self.types[name] = self
        self.efficiency = efficiency

    def __repr__(self):
        return "<{} name='{}'>".format(self.__class__.__name__, self.name)


MoleculeType(Constants.FRAGMENT, (False, False), True, 0.0, 0.98)
MoleculeType(Constants.PCR_PRODUCT, (True, True), False, 10.0, 0.95)
MoleculeType(Constants.PCR_PRODUCT_WITH_PRIMERS, (False, False), False, 10.0, 0.95)
MoleculeType(Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, (True, False), False, 10.0, 0.95)
MoleculeType(Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, (False, True), False, 10.0, 0.95)
MoleculeType(Constants.OVERLAP, (False, False), False, 0.0, 1.0)
MoleculeType(Constants.GAP, (False, False), False, 0.0, 1.0, synthesize=True)
