import bisect
import itertools
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
from Bio.SeqRecord import SeqRecord
from more_itertools import partition
from pyblast.utils import is_circular

from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.log import logger
from dasi.models import AlignmentContainer
from dasi.models import AlignmentGroup
from dasi.models import AssemblyNode
from dasi.models import MoleculeType
from dasi.models import PCRProductAlignmentGroup
from dasi.utils import bisect_between
from dasi.utils import Region
from dasi.utils import sort_with_keys
from dasi.utils.sequence_complexity import DNAStats


class AssemblyGraphBuilder:
    """Class that builds an AssemblyGraph from an alignment container."""

    COST_THRESHOLD = 10000

    def __init__(self, alignment_container: AlignmentContainer, span_cost=None):
        self.container = alignment_container
        if span_cost is None:
            self.span_cost = SpanCost.open()
        else:
            self.span_cost = span_cost
        self.G = None
        self.logger = logger(self)

    def add_node(self, node: AssemblyNode) -> None:
        """Add node to the graph.

        :param node: the assembly node to add
        :return: None
        """
        self.G.add_node(AssemblyNode(*node))

    def add_edge(
        self,
        n1: AssemblyNode,
        n2: AssemblyNode,
        cost: Union[float, None],
        material: Union[float, None],
        time: Union[float, None],
        efficiency: Union[float, None],
        span: int,
        atype: str,
        internal_or_external: str,
        condition: Tuple[bool, bool],
    ):
        """Add an edge between two assembly nodes.

        :param n1: src node
        :param n2: dest node
        :param name: name of the edge
        :param cost: overall cost of the edge
        :param material: material cost of the edge. Used in path calculations.
        :param time: time cost of the edge. Used in path calculations.
        :param efficiency: efficiency of the edge. Used in path calculations.
        :param span: spanning distance (in bp) of the edge.
        :param atype: alignment type of the edge.
        :param kwargs: additional kwargs for the edge data
        :return:
        """
        assert condition == atype.design
        assert internal_or_external == atype.int_or_ext
        self.G.add_edge(
            n1,
            n2,
            cost=cost,
            material=material,
            time=time,
            efficiency=efficiency,
            span=span,
            type_def=atype,
        )

    def iter_internal_edge_data(
        self, align: Union[AlignmentGroup, PCRProductAlignmentGroup]
    ) -> dict:
        q = align.query_region

        mtype = MoleculeType.types[align.type]
        a_expand, b_expand = mtype.design
        internal_cost = mtype.cost
        internal_efficiency = mtype.efficiency

        if q.cyclic:
            # cyclic
            if q.b < q.a:
                pairs = [(q.a, q.b + q.context_length), (q.a + q.context_length, q.b)]
            else:
                pairs = [(q.a, q.b), (q.a + q.context_length, q.b + q.context_length)]
        else:
            # linear
            pairs = [(q.a, q.b)]

        # nodes = []
        # edges = []
        for a, b in pairs:
            if a is not None:
                anode = (a, a_expand, "A")
                # nodes.append(anode)
                if b is not None:
                    bnode = (b, b_expand, "B")
                    # nodes.append(bnode)
                    yield (
                        anode,
                        bnode,
                        dict(
                            material=internal_cost,
                            cost=internal_cost / internal_efficiency,
                            time=0.1,
                            internal_or_external="internal",
                            span=len(align.query_region),
                            atype=MoleculeType.types[align.type],
                            efficiency=internal_efficiency,
                            condition=(a_expand, b_expand),
                        ),
                    )
                    # edges.append(edge)
        # return nodes, edges

    def add_internal_edges(
        self, groups: List[Union[AlignmentGroup, PCRProductAlignmentGroup]]
    ):
        for g in groups:
            for a, b, ab_data in self.iter_internal_edge_data(g):
                for a_overhang, b_overhang in itertools.product(
                    [True, False], repeat=2
                ):
                    a_node = AssemblyNode(a[0], a[1], a[2], a_overhang)
                    b_node = AssemblyNode(b[0], b[1], b[2], b_overhang)
                    self.add_edge(a_node, b_node, **ab_data)

    def add_external_edges(self, groups, group_keys, nodes: Iterable[AssemblyNode]):
        if not groups:
            return
        query_region = groups[0].query_region
        length = query_region.context_length
        a_nodes, b_nodes = partition(lambda x: x.type == "B", nodes)

        a_nodes = sorted(a_nodes, key=lambda x: x.index)
        b_nodes = sorted(b_nodes, key=lambda x: x.index)

        a_nodes_gap, a_nodes_overhang = [
            list(x) for x in partition(lambda x: x.overhang, a_nodes)
        ]
        b_nodes_gap, b_nodes_overhang = [
            list(x) for x in partition(lambda x: x.overhang, b_nodes)
        ]

        def make_overlap_iterator(a, b):
            """Find all nodes that satisfy the below condition: b.

            |--------|
                   |-------|
                   a

            With the exception that when a.index == b.index, this is not
            considered an overlap, but covered in the `make_gap_iterator`,
            due to using bisect.bisect_right. Overlap is more computationally
            intensive.
            """
            a, akeys = sort_with_keys(a, lambda x: x.index)
            b, bkeys = sort_with_keys(b, lambda x: x.index)
            for _a in a:
                i = bisect.bisect_right(bkeys, _a.index)
                for _b in b[i:]:
                    yield _b, _a

        def make_gap_itererator(a, b):
            """Find all nodes that satisfy the below condition: b.

            |--------|              |-------|              a
            """
            a, akeys = sort_with_keys(a, lambda x: x.index)
            b, bkeys = sort_with_keys(b, lambda x: x.index)

            for _a in a:
                i = bisect.bisect_right(bkeys, _a.index)
                for _b in b[:i]:
                    yield _b, _a

        def make_origin_iterator(a, b):
            a, akeys = sort_with_keys(a, lambda x: x.index)
            b, bkeys = sort_with_keys(b, lambda x: x.index)

            i = bisect.bisect_right(akeys, length)
            j = bisect.bisect_left(bkeys, length)
            for _a in a[:i]:
                for _b in b[j:]:
                    yield _b, _a

        overlap_iter = make_overlap_iterator(a_nodes_overhang, b_nodes_overhang)
        gap_iter = make_gap_itererator(a_nodes_gap, b_nodes_gap)
        gap_origin_iter = make_origin_iterator(a_nodes_gap, b_nodes_gap)
        overlap_origin_iter = make_origin_iterator(a_nodes_overhang, b_nodes_overhang)

        for bnode, anode in overlap_iter:
            self.add_overlap_edge(bnode, anode, query_region, group_keys, groups)

        for bnode, anode in gap_iter:
            self.add_gap_edge(bnode, anode, query_region)

        for bnode, anode in overlap_origin_iter:
            self.add_overlap_edge(
                bnode, anode, query_region, group_keys, groups, origin=True
            )
        for bnode, anode in gap_origin_iter:
            self.add_gap_edge(bnode, anode, query_region, origin=True)
        self._batch_add_edge_costs()

    def add_overlap_edge(
        self, bnode, anode, query_region, group_keys, groups, origin=False
    ):
        q = query_region.new(anode.index, bnode.index)
        if len(q) == 0:
            return
        if not origin:
            i, j = bisect_between(group_keys, q.a, q.b)
            filtered_groups = groups[i:j]
        else:
            i = bisect.bisect_left(group_keys, q.a)
            j = bisect.bisect_right(group_keys, q.b)
            filtered_groups = groups[i:] + groups[:j]
        if filtered_groups:
            if not origin:
                span = anode.index - bnode.index
            else:
                span = -len(q)
            if span <= 0:
                condition = (bnode.expandable, anode.expandable)
                self.add_edge(
                    bnode,
                    anode,
                    cost=None,
                    material=None,
                    time=None,
                    efficiency=None,
                    internal_or_external="external",
                    atype=MoleculeType.types[Constants.OVERLAP](condition),
                    condition=condition,
                    span=span,
                )

    def add_gap_edge(self, bnode, anode, query_region, origin=False):
        if not origin:
            span = anode.index - bnode.index
        else:
            span = len(query_region.new(bnode.index, anode.index))
        if span >= 0:
            condition = (bnode.expandable, anode.expandable)
            # cost_dict = self._get_cost(span, condition)
            self.add_edge(
                bnode,
                anode,
                cost=None,
                material=None,
                time=None,
                efficiency=None,
                internal_or_external="external",
                atype=MoleculeType.types[Constants.GAP](condition),
                span=span,
                condition=condition,
            )

    def _batch_add_edge_costs(self):
        """Add costs to all edges at once batch.

        :return:
        """
        # add external_edge_costs
        edge_dict = {}
        for n1, n2, edata in self.G.edges(data=True):
            if edata["cost"] is None:
                condition = edata["type_def"].design
                edge_dict.setdefault(condition, []).append(
                    ((n1, n2), edata, edata["span"])
                )
        edges_to_remove = []

        for condition, info in edge_dict.items():
            edges, edata, spans = zip(*info)
            npdf = self.span_cost.cost(np.array(spans), condition)
            data = npdf.aggregate(np.vstack)
            cost_i = [i for i, col in enumerate(npdf.columns) if col == "cost"][0]
            self.logger.debug(data.shape)
            # update each edge

            for i, (e, edge) in enumerate(zip(edata, edges)):
                if data[cost_i, i] > self.COST_THRESHOLD:
                    edges_to_remove.append(edge)
                else:
                    _d = {c: data[col_i, i] for col_i, c in enumerate(npdf.columns)}
                    e.update(_d)
        self.G.remove_edges_from(edges_to_remove)

    def _get_cost(self, bp, ext):
        data = self.span_cost(bp, ext).data
        return {k: v[0] for k, v in data.items()}

    def build_assembly_graph(self) -> nx.DiGraph:

        self.G = nx.DiGraph(name="Assembly Graph")

        groups, group_keys = sort_with_keys(
            self.container.get_groups_by_types(
                [
                    Constants.PCR_PRODUCT,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_PRIMERS,
                    Constants.PRIMER_EXTENSION_PRODUCT,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.FRAGMENT,
                ]
            ),
            key=lambda g: g.query_region.a,
        )

        self.add_internal_edges(groups)
        self.add_external_edges(groups, group_keys, self.G.nodes())

        nx.freeze(self.G)
        return self.G


class AssemblyGraphPostProcessor:
    def __init__(self, graph: nx.DiGraph, query: SeqRecord):
        self.graph = graph
        self.query = query
        self.stats = DNAStats(
            self.query + self.query,
            repeat_window=14,
            stats_window=20,
            hairpin_window=20,
        )
        self.logged_msgs = []
        self.COMPLEXITY_THRESHOLD = 12.0

    def update_edge_complexity(self, edata, complexity):
        edata["complexity"] = complexity
        if complexity > self.COMPLEXITY_THRESHOLD:
            edata["efficiency"] = 0.1
            return True
        return False

    # @staticmethod
    # def edge_to_region_helper(n1: int, n2: int, context_length: int, span: int, cyclic: bool):
    #     if span > 0:
    #         return Region(n1.index, n2.index, context_length, cyclic=cyclic)
    #     else:
    #         return Region(n2.index, n1.indes, context_length, cyclic=cyclic)
    #
    # def edge_to_region(self, n1, n2, span):
    #     cyclic = is_circular(self.query)
    #     length = len(self.query)
    #     return self.edge_to_region_helper(n1, n2, length, span, cyclic)

    # TODO: move this to a new class
    def complexity_update(self) -> nx.DiGraph:
        """Updates any gaps using the complexity measurements."""
        g = self.graph
        query = self.query
        bad_edges = []
        for n1, n2, edata in g.edges(data=True):
            if n1.type == "B" and n2.type == "A":
                span = edata["span"]
                if span > 0:
                    # TODO: cyclic may not always be tru
                    r = Region(
                        n1.index, n2.index, len(query), cyclic=is_circular(query)
                    )
                    score = self.stats.cost(r.a, r.c)
                    if self.update_edge_complexity(edata, score) is True:
                        bad_edges.append((n1, n2, edata))
                        self.logged_msgs.append("High complexity!")
        return bad_edges

    @staticmethod
    def optimize_partition(
        signatures: np.ndarray, step: int, i: int = None, j: int = None
    ):
        """Optimize partition by minimizing the number of signatures in the
        given array.

        :param signatures: array of signatures
        :param step: step size
        :param i:
        :param j:
        :return:
        """
        d = []

        if i is None:
            i = 0
        if j is None:
            j = signatures.shape[1]

        for x in range(i, j, step):
            m1 = np.empty(signatures.shape[1])
            m2 = m1.copy()
            m1.fill(np.nan)
            m2.fill(np.nan)

            m1[:x] = np.random.uniform(1, 10)
            m2[x:] = np.random.uniform(1, 10)

            d += [m1, m2]
        d = np.vstack(d)
        z = np.tile(d, signatures.shape[0]) * signatures.flatten()

        partition_index = np.repeat(
            np.arange(0, signatures.shape[1], step),
            signatures.shape[0] * signatures.shape[1] * 2,
        )

        a, b, c = np.unique(z, return_counts=True, return_index=True)
        i = b[np.where(c > 1)]
        a, c = np.unique(partition_index[i], return_counts=True)
        if len(c):
            arg = c.argmin()
            return a[arg], c[arg]

    # TODO: add logging to graph post processor
    # TODO: partition gaps
    def __call__(self):
        self.complexity_update()
        # if bad_edges:
        #     fwd = self.stats.fwd_signatures
        #     rev = self.stats.rev_signatures
        #     for n1, n2, edata in bad_edges:
        #         r = Region(
        #             n1.index, n2.index, len(self.query), cyclic=is_circular(self.query)
        #         )
        #         fwd_slice = fwd[r.a : r.c]
        #         rev_slice = rev[r.a : r.c]
        #         sig = np.vstack((fwd_slice, rev_slice))
        #         self.optimize_partition(sig, step=10)
