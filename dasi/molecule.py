from copy import copy
from typing import Tuple

from numpy import inf

from .constants import Constants


class MoleculeType:
    """Molecule metatype."""

    types = {}

    def __init__(
        self,
        name: str,
        design: Tuple[bool, bool],
        use_direct: bool,
        cost: float,
        efficiency=1.0,
        synthesize=False,
        min_size: int = None,
        max_size: int = None,
    ):
        self.name = name
        self.cost = cost
        self.use_direct = use_direct
        self.synthesize = synthesize
        self.types[name] = self
        self.efficiency = efficiency
        self.design = design
        self.int_or_ext = None
        self.min_size = min_size
        self.max_size = max_size

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
        min_size: int = None,
        max_size: int = None,
    ):
        super().__init__(
            name,
            use_direct=use_direct,
            design=design,
            cost=cost,
            efficiency=efficiency,
            synthesize=synthesize,
            min_size=min_size,
            max_size=max_size,
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
        min_size: int = None,
        max_size: int = None,
    ):
        super().__init__(
            name,
            use_direct=use_direct,
            design=None,
            cost=cost,
            efficiency=efficiency,
            synthesize=synthesize,
            min_size=min_size,
            max_size=max_size,
        )
        self.int_or_ext = "external"

    def __call__(self, design):
        copied = copy(self)
        copied.design = design
        return copied


InternalType(Constants.FRAGMENT, (False, False), True, 0.0, 0.98)
InternalType(Constants.PCR_PRODUCT, (True, True), False, 10.0, 0.95, min_size=100)
InternalType(
    Constants.PCR_PRODUCT_WITH_PRIMERS, (False, False), False, 10.0, 0.95, min_size=100
)
InternalType(
    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
    (True, False),
    False,
    10.0,
    0.95,
    min_size=100,
)
InternalType(
    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
    (False, True),
    False,
    10.0,
    0.95,
    min_size=100,
)
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

MoleculeType(
    name=Constants.TEMPLATE,
    design=None,
    use_direct=False,
    cost=0.0,
    efficiency=0.0,
    synthesize=False,
)


class Molecule:
    """An instance of a molecule type, with a sequence and which alignments are
    assigned to it."""

    def __init__(
        self, molecule_type, alignment_group, sequence, query_region=None, metadata=None
    ):
        self.type = molecule_type
        self.alignment_group = alignment_group
        self.sequence = sequence
        self.metadata = metadata or {}
        self.query_region = query_region

    def __repr__(self):
        return "<{cls} name='{name}' group='{group}'>".format(
            cls=self.__class__.__name__, name=self.type.name, group=self.alignment_group
        )


class Reaction:
    """An activity that takes in several Molecules and produces other
    Molecules."""

    def __init__(self, name, inputs, outputs):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs

    def __repr__(self):
        return "<{cls} name='{name}' outputs={products} regions={outputs}>".format(
            cls=self.__class__.__name__,
            name=self.name,
            products=[m.type.name for m in self.outputs],
            outputs=[m.query_region for m in self.outputs],
        )
