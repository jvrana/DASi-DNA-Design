from __future__ import annotations

from itertools import zip_longest
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from pyblast.utils import is_circular

from dasi.alignments import AlignmentContainer
from dasi.alignments import PCRProductAlignmentGroup
from dasi.constants import Constants
from dasi.constants import MoleculeType
from dasi.design.graph_builder import AssemblyNode
from dasi.exceptions import DasiDesignException
from dasi.log import logger
from dasi.utils import sort_cycle


class Assembly(Iterable):
    """Should take in a path, graph, container, seqdb to produce relevant
    information."""

    def __init__(
        self,
        nodes: List[AssemblyNode],
        container: AlignmentContainer,
        full_assembly_graph: nx.DiGraph,
        query_key: str,
        query: SeqRecord,
    ):
        self.logger = logger(self)
        self.container = container
        self.groups = container.groups()
        if len(self.groups) == 0:
            raise DasiDesignException("No groups were found in container.")
        self.query_key = query_key
        self.query = query
        self._nodes = tuple(nodes)
        self._full_graph = full_assembly_graph
        self.graph = self._subgraph(self._full_graph, nodes)
        nx.freeze(self.graph)

    @staticmethod
    def _missing_edata():
        return {
            "cost": np.inf,
            "weight": np.inf,
            "material": np.inf,
            "efficiency": 0.0,
            "type_def": MoleculeType.types[Constants.MISSING],
            "span": np.inf,
        }

    def _subgraph(self, graph: nx.DiGraph, nodes: List[AssemblyNode]):
        def _resolve(node: AssemblyNode, qregion) -> Tuple[AssemblyNode, dict]:
            new_node = AssemblyNode(qregion.t(node.index), *list(node)[1:])
            return new_node, {}

        subgraph = nx.OrderedDiGraph()
        nodes = [AssemblyNode(*n) for n in nodes]
        example_query_region = self.container.alignments[0].query_region

        resolved_nodes = [_resolve(node, example_query_region) for node in nodes]
        if self.cyclic:
            resolved_nodes = sort_cycle(
                resolved_nodes, key=lambda n: (n[0].type, n[0].index, n)
            )
        subgraph.add_nodes_from(resolved_nodes)

        if self.cyclic:
            pair_iter = list(pairwise(nodes + nodes[:1]))
        else:
            pair_iter = list(pairwise(nodes))

        for n1, n2 in pair_iter:
            edata = graph.get_edge_data(n1, n2)
            if edata is None:
                edata = self._missing_edata()
            else:
                assert edata["type_def"].int_or_ext

            # TODO: fix query_region (overlaps are backwards)
            query_region = self.container.alignments[0].query_region.new(
                n1.index, n2.index
            )
            groups = self.container.find_groups_by_pos(
                query_region.a,
                query_region.b,
                group_type=edata["type_def"].name,
                groups=self.groups,
            )
            if edata["type_def"].int_or_ext == "internal" and not groups:
                raise DasiDesignException(
                    "Missing groups for edge between {} and {}".format(n1, n2)
                )

            edata["groups"] = groups
            edata["query_region"] = query_region
            subgraph.add_edge(
                _resolve(n1, query_region)[0], _resolve(n2, query_region)[0], **edata
            )
        return subgraph

    @property
    def cyclic(self):
        return is_circular(self.query)

    # TODO: consolidate this with shortest path utils in networkx
    def cost(self):
        material = 0
        efficiency = 1.0
        for _, _, edata in self.edges():
            material += edata["material"]
            efficiency *= edata["efficiency"]
        if efficiency == 0:
            return np.inf
        return material / efficiency

    def edges(self, data=True) -> Iterable[Tuple[AssemblyNode, AssemblyNode, Dict]]:
        return self.graph.edges(data=data)

    def nodes(self, data=True) -> Iterable[Tuple[AssemblyNode, Dict]]:
        return self.graph.nodes(data=data)

    def edit_distance(
        self, other: Assembly, explain=False
    ) -> Union[int, Tuple[int, List[Tuple[int, str]]]]:
        differences = []
        for i, (n1, n2) in enumerate(
            zip_longest(self.nodes(data=False), other.nodes(data=False))
        ):
            if n1 is None or n2 is None:
                differences.append((i, "{} != {}".format(n1, n2)))
                continue
            if n1.index != n2.index:
                differences.append((i, "Index: {} != {}".format(n1.index, n2.index)))
            if n1.expandable != n2.expandable:
                differences.append(
                    (i, "Expandable: {} != {}".format(n1.expandable, n2.expandable))
                )
            if n1.type != n2.type:
                differences.append((i, "Type: {} != {}".format(n1.type, n2.type)))
            if n1.overhang != n2.overhang:
                differences.append(
                    (i, "Overhang: {} != {}".format(n1.overhang, n2.overhang))
                )
        dist = len(differences)
        if explain:
            return dist, differences
        return dist

    def print(self):
        print("query_name: {}".format(self.query.name))
        print("query_key: {}".format(self.query_key))
        print("Cost: {}".format(self.cost()))
        df = self.to_df()
        print(df)

    def print_diff(self, other: Assembly):
        for i, (n1, n2) in enumerate(
            zip_longest(self.nodes(data=False), other.nodes(data=False))
        ):
            if n1 != n2:
                desc = False
            else:
                desc = True
            print("{} {} {}".format(desc, n1, n2))

    def to_df(self):
        rows = []
        for n1, n2, edata in self.edges():
            groups = edata["groups"]

            if groups:
                group = groups[0]
                if isinstance(group, PCRProductAlignmentGroup):
                    alignments = group.alignments
                else:
                    alignments = group.alignments[:1]
            else:
                alignments = []
            subject_keys = [a.subject_key for a in alignments]
            subject_names = [self.container.seqdb[key].name for key in subject_keys]
            subject_starts = [a.subject_region.a for a in alignments]
            subject_ends = [a.subject_region.b for a in alignments]

            rows.append(
                {
                    "query_start": edata["query_region"].a,
                    "query_end": edata["query_region"].b,
                    "subject_names": subject_names,
                    "subject_keys": subject_keys,
                    "subject_start": subject_starts,
                    "subject_ends": subject_ends,
                    "cost": edata["cost"],
                    "material": edata["material"],
                    "span": edata["span"],
                    "type": edata["type_def"].name,
                    "internal_or_external": edata["type_def"].int_or_ext,
                    "efficiency": edata.get("efficiency", np.nan),
                }
            )

        df = pd.DataFrame(rows)
        return df

    def __eq__(self, other: Assembly) -> bool:
        return self.edit_distance(other) == 0

    def __iter__(self) -> Generator[AssemblyNode]:
        for n in self.nodes(data=False):
            yield n