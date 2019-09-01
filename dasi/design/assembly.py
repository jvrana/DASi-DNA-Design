from dasi.cost import SpanCost
from dasi.log import logger
from dasi.constants import Constants
from dasi.alignments import AlignmentContainer
from dasi.utils import sort_with_keys, bisect_slice_between
import itertools
import networkx as nx
from collections import namedtuple

from more_itertools import partition, unique_everseen

AssemblyNode = namedtuple('AssemblyNode', 'index expandable type overhang')


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

    def add_edge(self, n1, n2, weight, name, span, type):
        self.G.add_edge(AssemblyNode(*n1), AssemblyNode(*n2), weight=weight, name=name, span=span, type=type)

    def internal_edge_data(self, align):
        q = align.query_region
        a_expand, b_expand = True, True
        if align.type in [Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            b_expand = False
        if align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            a_expand = False

        # TODO: internal edge cost?
        ### INTERNAL EDGE
        if align.type == Constants.FRAGMENT:
            internal_cost = 0
        elif align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER]:
            internal_cost = 30 + 30
        elif align.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
            internal_cost = 30
        elif align.type == Constants.PCR_PRODUCT:
            internal_cost = 30 + 30 + 30

        a = q.a
        if q.b < q.a and q.cyclic:
            # is cyclic
            b = q.b + q.context_length
            lengths = [0]
        else:
            # is linear
            lengths = [0, q.context_length]
            b = q.b

        nodes = []
        for length in lengths:
            node = (a + length, a_expand, 'A'), (b + length, b_expand, 'B'), dict(
                            weight=internal_cost,
                            name='',
                            span=len(align.query_region),
                            type=align.type
            )
            nodes.append(node)
        return nodes

    def add_internal_edges(self, groups):
        for g in groups:
            for a, b, ab_data in self.internal_edge_data(g):
                for a_overhang, b_overhang in itertools.product([True, False], repeat=2):
                    a_node = (a[0], a[1], a[2], a_overhang)
                    b_node = (b[0], b[1], b[2], b_overhang)
                    self.add_edge(a_node, b_node, **ab_data)

    def add_external_edges(self, groups, group_keys):
        if not groups:
            return
        query_region = groups[0].query_region
        a_nodes, b_nodes = partition(lambda x: x.type == 'B', self.G.nodes())
        for (b, b_expand, bid, b_overhang), (a, a_expand, aid, a_overhang) in itertools.product(b_nodes, a_nodes):
            if not (b_overhang or a_overhang):
                ba = query_region.new(b, a)

                # TODO: PRIORITY no way to determine overlaps from just end points
                cost, desc = self.span_cost.cost_and_desc(len(ba), (b_expand, a_expand))
                if cost < self.COST_THRESHOLD:
                    n1 = (b, b_expand, bid, b_overhang)
                    n2 = (a, a_expand, aid, a_overhang)
                    self.add_edge(
                        n1, n2,
                        weight=cost,
                        name='',
                        span=len(ba),
                        type=desc
                    )
            else:
                # TODO: PRIORITY this step is extremely slow
                filtered_groups = bisect_slice_between(groups, group_keys, a, b)
                if filtered_groups:
                    ab = filtered_groups[0].query_region.new(a, b)
                    cost, desc = self.span_cost.cost_and_desc(-len(ab), (b_expand, a_expand))

                    n1 = (b, b_expand, bid, True)
                    n2 = (a, a_expand, aid, True)

                    if cost < self.COST_THRESHOLD:
                        self.add_edge(n1, n2, weight=cost, name='overlap', span=-len(ab), type=desc)

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
        self.add_external_edges(groups, group_keys)
        return self.G
