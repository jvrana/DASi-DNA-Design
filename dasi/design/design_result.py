import bisect
from typing import Generator
from typing import Iterable
from typing import List
from typing import Tuple

import networkx as nx

from dasi.exceptions import DasiInvalidMolecularAssembly
from dasi.models import AlignmentContainer
from dasi.models import Assembly
from dasi.models import AssemblyNode


class DesignResult(Iterable):
    """DesignResult container.

    Maintains a list of top assemblies.
    """

    def __init__(
        self, container: AlignmentContainer, graph: nx.DiGraph, query_key: str
    ):
        self.container = container
        self.graph = graph
        self.query_key = query_key
        self.query = self.container.seqdb[query_key]
        self._assemblies = []
        self._keys = []

    @property
    def seqdb(self):
        return self.container.seqdb

    @property
    def assemblies(self) -> Tuple[Assembly, ...]:
        """Return a tuple of all assemblies.

        :return: tuple of all assemblies.
        """
        return tuple(self._assemblies)

    def _add_assembly_from_path(self, path: List[AssemblyNode]):
        return Assembly(
            path,
            self.container,
            self.graph,
            self.query_key,
            self.query,
            seqdb=self.seqdb,
            do_raise=False,
        )

    def add_assembly(
        self,
        path: List[AssemblyNode],
        ignore_invalid: bool = False,
        allow_invalid: bool = False,
    ):
        """Add an assembly from a list of nodes.

        :param path: list of nodes
        :return: None
        """

        assembly = self._add_assembly_from_path(path)

        # Validate the assembly
        try:
            assembly.post_validate()
        except DasiInvalidMolecularAssembly as e:
            if ignore_invalid:
                return
            elif not allow_invalid:
                raise e

        cost = assembly.compute_cost()
        n_nodes = len(assembly._nodes)
        k = (cost, n_nodes)
        i = bisect.bisect_left(self._keys, k)
        self._assemblies.insert(i, assembly)
        self._keys.insert(i, k)
        return assembly

    def add_assemblies(
        self,
        paths: List[List[AssemblyNode]],
        ignore_invalid: bool = False,
        allow_invalid: bool = False,
    ):
        """Adds a list of assemblies.

        :param paths: list of list of paths
        :return: None
        """
        for path in paths:
            self.add_assembly(
                path, ignore_invalid=ignore_invalid, allow_invalid=allow_invalid
            )
    #
    #
    #     for a in self.assemblies:
    #         for n1, n2, edata in a.edges():
    #             result = edata["sequence_result"]
    #             mol_type = edata["type_def"]
    #             if mol_type.int_or_ext == "internal":
    #                 pass
    #                 # add a new reaction
    #                 # add template
    #                 # add primers
    #             elif mol_type.use_direct:
    #                 pass
    #                 # add fragment
    #             elif mol_type.synthesize:
    #                 pass
    #                 # add synthesize
    #             else:
    #                 pass
    #                 # raise Exception

    def __iter__(self) -> Generator[Assembly, None, None]:
        """Yield assemblies.

        :yield: assembly
        """
        yield from self.assemblies

    def __getitem__(self, item: str) -> Assembly:
        return list(self)[item]

    def __str__(self):
        return "<{cls}, query={qname} {qk} nassemblies={n}>".format(
            cls=self.__class__.__name__,
            qname=self.query.name,
            qk=self.query_key,
            n=len(self.assemblies),
        )
