from dasi.cost import SpanCost
from dasi.log import logger
from dasi.constants import Constants
from dasi.alignments import AlignmentContainer
from dasi.utils import sort_with_keys, bisect_slice_between
from bisect import bisect_left
import itertools
import networkx as nx


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

    def alignment_to_internal_edge(self, align):
        q = align.query_region
        a_expand, b_expand = True, True
        if align.type in [Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            b_expand = False
        if align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.FRAGMENT, Constants.PCR_PRODUCT_WITH_PRIMERS]:
            a_expand = False

        ### INTERNAL EDGE
        if align.type == Constants.FRAGMENT:
            internal_cost = 0
        elif align.type in [Constants.PCR_PRODUCT_WITH_LEFT_PRIMER, Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER]:
            internal_cost = 30 + 30
        elif align.type == Constants.PCR_PRODUCT_WITH_PRIMERS:
            internal_cost = 30
        elif align.type == Constants.PCR_PRODUCT:
            internal_cost = 30 + 30 + 30

        return (q.a, a_expand, 'A'), (q.b, b_expand, 'B'), dict(
                    weight=internal_cost,
                    name='',
                    span_length=len(align.query_region),
                    type=align.type
                )

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

        # TODO: reduce number of times edge

        a_arr = set()
        b_arr = set()

        for g in groups:
            a, b, ab_data = self.alignment_to_internal_edge(g)
            for a_overhang, b_overhang in itertools.product([True, False], repeat=2):
                a_node = (a[0], a[1], a[2], a_overhang)
                b_node = (b[0], b[1], b[2], b_overhang)
                self.G.add_edge(a_node, b_node, **ab_data)
                a_arr.add(a_node)
                b_arr.add(b_node)

        ### EXTERNAL EDGES
        if groups:
            query = groups[0].query_region
            for (b, b_expand, bid, b_overhang), (a, a_expand, aid, a_overhang) in itertools.product(b_arr, a_arr):
                if not (b_overhang or a_overhang):
                    ba = query.new(b, a)
                    # ab = query.new(a, b)

                    # TODO: PRIORITY no way to determine overlaps from just end points
                    cost, desc = self.span_cost.cost_and_desc(len(ba), (b_expand, a_expand))
                    if cost < self.COST_THRESHOLD:
                        self.G.add_edge(
                            (b, b_expand, bid, b_overhang),
                            (a, a_expand, aid, a_overhang),
                            weight=cost,
                            name='',
                            span_length=len(ba),
                            type=desc
                        )
                else:
                    pass
                    # # TODO: PRIORITY this step is extremely slow
                    # filtered_groups = bisect_slice_between(groups, group_keys, a, b)
                    # if filtered_groups:
                    #     ab = filtered_groups[0].query_region.new(a, b)
                    #     cost, desc = self.span_cost.cost_and_desc(-len(ab), (b_expand, a_expand))
                    #     if cost < self.COST_THRESHOLD:
                    #             self.G.add_edge(
                    #                 (b, b_expand, bid, True),
                    #                 (a, a_expand, aid, True),
                    #                 weight=cost,
                    #                 name='overlap',
                    #                 span_length=-len(ab),
                    #                 type=desc
                    #         )
        return self.G
