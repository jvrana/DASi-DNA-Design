from dasi.cost import SpanCost
from dasi.log import logger
from dasi.constants import Constants
from dasi.alignments import AlignmentContainer
from dasi.utils import sort_with_keys, bisect_slice_between
import itertools
import networkx as nx
from collections import namedtuple
from dasi.exceptions import DASiException
from more_itertools import partition
from typing import Iterable
import bisect

AssemblyNode = namedtuple("AssemblyNode", "index expandable type overhang")


class AssemblyGraphBuilder(object):
    """
    Class that builds an AssemblyGraph from an alignment container.
    """

    COST_THRESHOLD = 10000

    def __init__(self, alignment_container: AlignmentContainer, span_cost=None):
        self.container = alignment_container
        if span_cost is None:
            self.span_cost = SpanCost()
        else:
            self.span_cost = span_cost
        self.G = None
        self.logger = logger(self)

    def add_node(self, node):
        self.G.add_node(AssemblyNode(*node))

    def add_edge(self, n1, n2, weight, name, span, type):
        n1 = AssemblyNode(*n1)
        n2 = AssemblyNode(*n2)
        self.G.add_edge(n1, n2, weight=weight, name=name, span=span, type=type)

    # TODO: internal cost should correlate with span_cost
    # TODO: interanl cost should also have an efficiency associated with it (can you amplify it?).
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

    @staticmethod
    def internal_edge_cost(align):
        if align.type == Constants.FRAGMENT:
            internal_cost = 0
        elif align.type in [
            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
            Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
        ]:
            internal_cost = 30 + 30
        elif align.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
            internal_cost = 30
        elif align.type == Constants.PCR_PRODUCT:
            internal_cost = 30 + 30 + 30
        else:
            raise DASiException("Could not determine cost of {}".format(align))
        return internal_cost

    def internal_edge_data(self, align):
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

        nodes = []
        edges = []
        for a, b in pairs:
            if a is not None:
                anode = (a, a_expand, "A")
                nodes.append(anode)
                if b is not None:
                    bnode = (b, b_expand, "B")
                    nodes.append(bnode)
                    edge = (
                        anode,
                        bnode,
                        dict(
                            weight=internal_cost,
                            name="",
                            span=len(align.query_region),
                            type=align.type,
                        ),
                    )
                    edges.append(edge)
        return nodes, edges

    def add_internal_edges(self, groups):
        for g in groups:
            nodes, edges = self.internal_edge_data(g)

            for a_overhang, b_overhang in itertools.product([True, False], repeat=2):
                # for n in nodes:
                #     if n[2] == 'A':
                #         o = a_overhang
                #     else:
                #         o = b_overhang
                #     self.add_node((n[0], n[1], n[2], o))
                for a, b, ab_data in edges:
                    a_node = (a[0], a[1], a[2], a_overhang)
                    b_node = (b[0], b[1], b[2], b_overhang)
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

        def make_overlap_iterator(anodes, bnodes):
            anodes, akeys = sort_with_keys(anodes, lambda x: x.index)
            for bnode in bnodes:
                i = bisect.bisect_left(akeys, bnode.index)
                for anode in anodes[:i]:
                    yield bnode, anode

        def make_gap_itererator(anodes, bnodes):
            anodes, akeys = sort_with_keys(anodes, lambda x: x.index)
            for bnode in bnodes:
                i = bisect.bisect_left(akeys, bnode.index)
                for anode in anodes[i:]:
                    yield bnode, anode

        def make_origin_iterator(anodes, bnodes):
            anodes, akeys = sort_with_keys(anodes, lambda x: x.index)
            bnodes, bkeys = sort_with_keys(bnodes, lambda x: x.index)

            i = bisect.bisect_right(akeys, length)
            j = bisect.bisect_left(bkeys, length)

            for bnode in bnodes[j:]:
                for anode in anodes[:i]:
                    yield bnode, anode

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

    def add_overlap_edge(
        self, bnode, anode, query_region, group_keys, groups, origin=False
    ):
        # TODO: PRIORITY this step is extremely slow
        q = query_region.new(anode.index, bnode.index, allow_wrap=True)
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
                cost, desc = self.span_cost.cost_and_desc(
                    span, (bnode.expandable, anode.expandable)
                )
                if cost < self.COST_THRESHOLD:
                    self.add_edge(
                        bnode, anode, weight=cost, name="overlap", span=span, type=desc
                    )

    def add_gap_edge(self, bnode, anode, query_region, origin=False):
        # TODO: PRIORITY no way to determine overlaps from just end points
        if not origin:
            span = anode.index - bnode.index
        else:
            span = len(query_region.new(bnode.index, anode.index))
        if span >= 0:
            cost, desc = self.span_cost.cost_and_desc(
                span, (bnode.expandable, anode.expandable)
            )
            if cost < self.COST_THRESHOLD:
                self.add_edge(bnode, anode, weight=cost, name="", span=span, type=desc)

    def build_assembly_graph(self):

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
