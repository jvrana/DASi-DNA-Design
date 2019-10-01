import bisect
import itertools
from collections import namedtuple
from typing import Iterable, Union, List

import networkx as nx
import numpy as np
from more_itertools import partition

from dasi.alignments import AlignmentContainer, AlignmentGroup, ComplexAlignmentGroup
from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.exceptions import DASiException
from dasi.log import logger
from dasi.utils import sort_with_keys, bisect_slice_between

AssemblyNode = namedtuple("AssemblyNode", "index expandable type overhang")


class AssemblyGraphBuilder(object):
    """
    Class that builds an AssemblyGraph from an alignment container.
    """

    COST_THRESHOLD = 10000

    def __init__(self, alignment_container: AlignmentContainer, span_cost=None):
        self.container = alignment_container
        if span_cost is None:
            self.span_cost = SpanCost.default()
        else:
            self.span_cost = span_cost
        self.G = None
        self.logger = logger(self)

    def add_node(self, node: AssemblyNode) -> None:
        """
        Add node to the graph.

        :param node: the assembly node to add
        :return: None
        """
        self.G.add_node(AssemblyNode(*node))

    def add_edge(
        self,
        n1: AssemblyNode,
        n2: AssemblyNode,
        name,
        cost: Union[float, None],
        material: Union[float, None],
        time: Union[float, None],
        efficiency: Union[float, None],
        span: int,
        atype: str,
        **kwargs
    ):
        """
        Add an edge between two assembly nodes.

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
        self.G.add_edge(
            n1,
            n2,
            name=name,
            cost=cost,
            material=material,
            time=time,
            efficiency=efficiency,
            span=span,
            type=atype,
            **kwargs
        )

    # TODO: internal cost should correlate with span_cost
    # TODO: internal cost should also have an efficiency associated with it (can you amplify it?).
    #       How fast can you compute this for all alignments?
    # TODO: test internal and external edge costs for all cases.
    # TODO: We are counting the primer material cost **twice**, once in the JxnCost and once here.
    #       Really, the cost of the internal fragment is either its already available ($0) or
    #       we have to PCR amplify. The material cost of PCR amplification is just the reagents and time
    #       of setting up the reaction itself.
    #       We should also account for the 'efficiency' of PCR amplification and include that in the edge
    #       Perhaps multiple edges and taking the 'max' efficiency and 'min' material cost, recording
    #       which alignment is being used for the calculations.
    #       For the jxn cost, we need a way to optimize the JxnEfficiency by *efficiency* AND *material cost*
    #       while heavily emphasizing the 'efficiency', since that has the most impact on overall assembly cost.
    #       At the end, we will need to separate material and efficiency and use efficiency in the overall assembly
    #       optimization.

    # TODO: make internal_cost, and efficiency a global parameter
    @staticmethod
    def internal_edge_cost(
        align: Union[AlignmentGroup, ComplexAlignmentGroup]
    ) -> float:
        """

        :param align:
        :return:
        """
        if align.type == Constants.FRAGMENT:
            internal_cost = 0
        elif align.type in [
            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
            Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
        ]:
            internal_cost = 10.0
        elif align.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
            internal_cost = 10.0
        elif align.type == Constants.PCR_PRODUCT:
            internal_cost = 10.0
        else:
            raise DASiException("Could not determine cost of {}".format(align))
        return internal_cost

    def iter_internal_edge_data(
        self, align: Union[AlignmentGroup, ComplexAlignmentGroup]
    ) -> dict:
        q = align.query_region
        a_expand, b_expand = True, True
        if align.type in [
            Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
            Constants.FRAGMENT,
            Constants.PCR_PRODUCT_WITH_PRIMERS,
        ]:
            b_expand = False
        if align.type in [
            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
            Constants.FRAGMENT,
            Constants.PCR_PRODUCT_WITH_PRIMERS,
        ]:
            a_expand = False

        internal_cost = self.internal_edge_cost(align)

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
                            weight=internal_cost / 0.95,
                            material=internal_cost,
                            cost=internal_cost / 0.95,
                            name="",
                            time=0.1,
                            internal_or_external="internal",
                            span=len(align.query_region),
                            condition=(a_expand, b_expand),
                            atype=align.type,
                            efficiency=0.95,
                        ),
                    )
                    # edges.append(edge)
        # return nodes, edges

    def add_internal_edges(
        self, groups: List[Union[AlignmentGroup, ComplexAlignmentGroup]]
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
            a, akeys = sort_with_keys(a, lambda x: x.index)
            for _b in b:
                i = bisect.bisect_left(akeys, _b.index)
                for _a in a[:i]:
                    yield _b, _a

        def make_gap_itererator(a, b):
            a, akeys = sort_with_keys(a, lambda x: x.index)
            for _b in b:
                i = bisect.bisect_left(akeys, _b.index)
                for _a in a[i:]:
                    yield _b, _a

        def make_origin_iterator(a, b):
            a, akeys = sort_with_keys(a, lambda x: x.index)
            b, bkeys = sort_with_keys(b, lambda x: x.index)

            i = bisect.bisect_right(akeys, length)
            j = bisect.bisect_left(bkeys, length)

            for _b in b[j:]:
                for _a in a[:i]:
                    yield _b, _a

        overlap_iter = make_overlap_iterator(a_nodes_overhang, b_nodes_overhang)
        gap_iter = make_gap_itererator(a_nodes_gap, b_nodes_gap)
        gap_origin_iter = make_origin_iterator(a_nodes_gap, b_nodes_gap)
        overlap_origin_iter = make_origin_iterator(a_nodes_overhang, b_nodes_overhang)

        for bnode, anode in overlap_iter:
            self.add_overlap_edge(bnode, anode, query_region, group_keys, groups)
        for bnode, anode in gap_iter:
            self.add_gap_edge(bnode, anode, query_region)

        for bnode, anode in gap_origin_iter:
            self.add_gap_edge(bnode, anode, query_region, origin=True)
        for bnode, anode in overlap_origin_iter:
            self.add_overlap_edge(
                bnode, anode, query_region, group_keys, groups, origin=True
            )
        self._batch_add_edge_costs()

    def _batch_add_edge_costs(self):
        # add external_edge_costs
        edge_dict = {}
        for n1, n2, edata in self.G.edges(data=True):
            if edata["cost"] is None:
                condition = edata["condition"]
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

    def add_overlap_edge(
        self, bnode, anode, query_region, group_keys, groups, origin=False
    ):
        # TODO: PRIORITY this step is extremely slow
        q = query_region.new(anode.index, bnode.index)
        if len(q) == 0:
            return
        if not origin:
            filtered_groups = bisect_slice_between(groups, group_keys, q.a, q.b)
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
                    # weight=cost_dict["cost"],
                    weight=None,
                    cost=None,
                    material=None,
                    time=None,
                    efficiency=None,
                    internal_or_external="external",
                    name="overlap",
                    atype="overlap",
                    condition=condition,
                    span=span,
                )

    def add_gap_edge(self, bnode, anode, query_region, origin=False):
        # TODO: PRIORITY no way to determine overlaps from just end points
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
                weight=None,
                cost=None,
                material=None,
                time=None,
                efficiency=None,
                internal_or_external="external",
                name="gap",
                atype="gap",
                condition=condition,
                span=span,
            )

    def build_assembly_graph(self) -> nx.DiGraph:

        self.G = nx.DiGraph(name="Assembly Graph")

        groups, group_keys = sort_with_keys(
            self.container.get_groups_by_types(
                [
                    Constants.PCR_PRODUCT,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_PRIMERS,
                    Constants.FRAGMENT,
                ]
            ),
            key=lambda g: g.query_region.a,
        )

        self.add_internal_edges(groups)
        self.add_external_edges(groups, group_keys, self.G.nodes())

        nx.freeze(self.G)
        return self.G
