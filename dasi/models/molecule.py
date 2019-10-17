"""classes representing *molecules* and types of *molecules*"""
from __future__ import annotations

from copy import copy
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import inf

from dasi.constants import Constants
from dasi.utils import Region


class MoleculeType:
    """Molecule metatype."""

    types = {}

    def __init__(
        self,
        name: str,
        design: Union[None, Tuple[bool, bool]],
        use_direct: bool,
        cost: float,
        efficiency=1.0,
        synthesize=False,
        min_size: int = None,
        max_size: int = None,
    ):
        """Initializes a new molecule."""
        self.name = name
        self.cost = cost
        self.use_direct = use_direct
        self.synthesize = synthesize
        if name in self.types:
            raise ValueError("Cannot re-define molecule type '{}'".format(name))
        self.types[name] = self
        self.efficiency = efficiency
        self.design = design
        self.int_or_ext = None
        self.min_size = min_size
        self.max_size = max_size

    def __repr__(self):
        return "<{} name='{}'>".format(self.__class__.__name__, self.name)


class InternalType(MoleculeType):
    """Molecule representing physical molecule that can be provided by the
    lab."""

    def __init__(
        self,
        name: str,
        design: Tuple[bool, bool],
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
    """Molecule to be designed."""

    def __init__(
        self,
        name: str,
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


# TODO: efficiencies and base costs could be parameters
# TODO: min_size and max_size could be parameters
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
        self,
        molecule_type: MoleculeType,
        alignment_group,
        sequence: SeqRecord,
        query_region: Region = None,
        metadata: Dict = None,
    ):
        self.type = molecule_type
        self.alignment_group = alignment_group
        assert issubclass(type(sequence), SeqRecord)
        assert issubclass(type(sequence.seq), Seq)
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

    def __init__(self, name: str, inputs: List[Molecule], outputs: List[Molecule]):
        self.name = name
        self.inputs = inputs  #: input molecules to the reaction
        self.outputs = outputs  #: output molecule of the reaction

    def __repr__(self):
        return "<{cls} name='{name}' outputs={products} regions={outputs}>".format(
            cls=self.__class__.__name__,
            name=self.name,
            products=[m.type.name for m in self.outputs],
            outputs=[m.query_region for m in self.outputs],
        )
