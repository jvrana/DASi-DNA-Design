import bisect
import itertools
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

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

    @staticmethod
    def add_edge_to_graph(
        graph,
        n1: AssemblyNode,
        n2: AssemblyNode,
        cost: Union[float, None],
        material: Union[float, None],
        time: Union[float, None],
        efficiency: Union[float, None],
        span: int,
        atype: MoleculeType,
        group: Union[AlignmentGroup, PCRProductAlignmentGroup],
        notes: str = "",
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
        return graph.add_edge(
            n1,
            n2,
            cost=cost,
            material=material,
            time=time,
            efficiency=efficiency,
            span=span,
            type_def=atype,
            group=group,
            notes=notes,
        )

    def add_edge(
        self,
        n1: AssemblyNode,
        n2: AssemblyNode,
        cost: Union[float, None],
        material: Union[float, None],
        time: Union[float, None],
        efficiency: Union[float, None],
        span: int,
        atype: MoleculeType,
        internal_or_external: str,
        condition: Tuple[bool, bool],
        group: Union[AlignmentGroup, PCRProductAlignmentGroup],
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
        return self.add_edge_to_graph(
            graph=self.G,
            n1=n1,
            n2=n2,
            cost=cost,
            material=material,
            time=time,
            efficiency=efficiency,
            span=span,
            atype=atype,
            group=group,
            notes="",
        )

    def iter_internal_edge_data(
        self, group: Union[AlignmentGroup, PCRProductAlignmentGroup]
    ) -> dict:
        q = group.query_region

        mtype = MoleculeType.types[group.type]
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
                            span=len(group.query_region),
                            atype=MoleculeType.types[group.type],
                            efficiency=internal_efficiency,
                            condition=(a_expand, b_expand),
                            group=group,
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
                    group=None,
                )

    def add_gap_edge(self, bnode, anode, query_region, origin=False):
        try:
            if not origin:
                span = anode.index - bnode.index
            else:
                if bnode.index == query_region.context_length:
                    span = len(query_region.new(0, anode.index))
                else:
                    span = len(query_region.new(bnode.index, anode.index))
        except Exception as e:
            print(
                "An exception was raised while attempting to create a abstract region"
            )
            print("origin == " + str(origin))
            print("bnode == " + str(bnode))
            print("anode == " + str(anode))
            print("query_region == " + str(query_region))
            raise e
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
                group=None,
            )

    @staticmethod
    def batch_add_edge_costs(graph, edges, span_cost, cost_threshold: float = None):
        """Add costs to all edges at once.

        :return:
        """
        # add external_edge_costs
        edge_dict = {}
        for n1, n2, edata in edges:
            condition = edata["type_def"].design
            assert isinstance(condition, tuple)
            edge_dict.setdefault(condition, []).append(((n1, n2), edata, edata["span"]))

        edges_to_remove = []
        for condition, info in edge_dict.items():
            edges, edata, spans = zip(*info)
            npdf = span_cost.cost(np.array(spans), condition)
            data = npdf.aggregate(np.vstack)
            cost_i = [i for i, col in enumerate(npdf.columns) if col == "cost"][0]

            for i, (e, edge) in enumerate(zip(edata, edges)):
                if cost_threshold is not None and data[cost_i, i] > cost_threshold:
                    edges_to_remove.append(edge)
                else:
                    _d = {c: data[col_i, i] for col_i, c in enumerate(npdf.columns)}
                    e.update(_d)
        if cost_threshold is not None:
            graph.remove_edges_from(edges_to_remove)

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
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_PRIMERS,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.FRAGMENT,
                    Constants.SHARED_SYNTHESIZED_FRAGMENT,
                ]
            ),
            key=lambda g: g.query_region.a,
        )

        self.add_internal_edges(groups)
        self.add_external_edges(groups, group_keys, self.G.nodes())

        edges = []
        for n1, n2, edata in self.G.edges(data=True):
            if edata["cost"] is None:
                edges.append((n1, n2, edata))
            elif edata["type_def"].name == Constants.SHARED_SYNTHESIZED_FRAGMENT:
                edges.append((n1, n2, edata))
        self.batch_add_edge_costs(self.G, edges, self.span_cost, self.COST_THRESHOLD)

        # TODO: freeze?
        # nx.freeze(self.G)
        return self.G


class AssemblyGraphPostProcessor:
    def __init__(self, graph: nx.DiGraph, query: SeqRecord, span_cost: SpanCost):
        self.graph = graph
        self.query = query
        self.stats = DNAStats(
            self.query + self.query,
            repeat_window=14,
            stats_window=20,
            hairpin_window=20,
        )
        self.logged_msgs = []
        self.COMPLEXITY_THRESHOLD = 10.0
        self.logger = logger(self)
        self.span_cost = span_cost

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

    def _edge_to_region(self, n1, n2):
        if n1.index == len(self.query) and is_circular(self.query):
            a = 0
            b = n2.index
        else:
            a = n1.index
            b = n2.index
        return Region(a, b, len(self.query), cyclic=is_circular(self.query))

    def _complexity_to_efficiency(self, edata):
        if edata["complexity"] > self.COMPLEXITY_THRESHOLD:
            edata["efficiency"] = 0.1
            return True
        return False

    def update_complexity(self, n1, n2, edata):
        bad_edges = []
        if edata["type_def"].synthesize:
            span = edata["span"]
            if span > 0:
                # TODO: cyclic may not always be tru
                r = self._edge_to_region(n1, n2)
                complexity = self.stats.cost(r.a, r.c)
                edata["complexity"] = complexity
                if self._complexity_to_efficiency(edata):
                    bad_edges.append((n1, n2, edata))
        return bad_edges

    def update_long_pcr_products(self, n1, n2, edata):
        if edata["type_def"].name in [
            Constants.PCR_PRODUCT,
            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
            Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
            Constants.PCR_PRODUCT_WITH_PRIMERS,
        ]:
            span = edata["span"]
            if span > 4000:
                edata["efficiency"] *= 0.8
            elif span > 5000:
                edata["efficiency"] *= 0.5

    def synthesis_partitioner(self, n1, n2, edata, border):
        r = self._edge_to_region(n1, n2)
        partitions = self.stats.partition(
            10,
            25,
            i=r.a,
            j=r.c,
            border=border,
            stopping_threshold=self.COMPLEXITY_THRESHOLD,
        )

        if not partitions:
            return []

        best_partition = partitions[0]
        n3 = AssemblyNode(
            best_partition["index_1"][1], False, str(uuid4()), overhang=True
        )
        n4 = AssemblyNode(
            best_partition["index_2"][0], False, str(uuid4()), overhang=True
        )

        AssemblyGraphBuilder.add_edge_to_graph(
            self.graph,
            n1,
            n3,
            cost=None,
            material=None,
            time=None,
            efficiency=None,
            span=len(self._edge_to_region(n1, n3)),
            atype=edata["type_def"],
            group=edata["group"],
        )

        # internal edge
        AssemblyGraphBuilder.add_edge_to_graph(
            self.graph,
            n3,
            n4,
            cost=None,
            material=None,
            time=None,
            efficiency=None,
            span=-len(self._edge_to_region(n4, n3)),
            atype=MoleculeType.types[Constants.OVERLAP],
            group=None,
        )

        # overlap edge
        AssemblyGraphBuilder.add_edge_to_graph(
            self.graph,
            n4,
            n2,
            cost=None,
            material=None,
            time=None,
            efficiency=None,
            span=-len(self._edge_to_region(n4, n2)),
            atype=edata["type_def"],
            group=edata["group"],
        )

        edata1 = self.graph.get_edge_data(n1, n3)
        edata2 = self.graph.get_edge_data(n3, n4)
        edata3 = self.graph.get_edge_data(n4, n2)

        edata1["complexity"] = best_partition["cost_1"]
        edata3["complexity"] = best_partition["cost_2"]
        edata3["notes"] += " Partitioned"
        edata2["notes"] += " Partitioned"
        edata1["notes"] += " Partitioned"

        edge1 = (n1, n3, edata1)
        edge2 = (n3, n4, edata2)
        edge3 = (n4, n2, edata3)
        return [edge1, edge2, edge3]

    def partition(self, edges):
        new_edges = []
        for n1, n2, edata in self.logger.tqdm(
            edges, "INFO", desc="Paritioning sequences..."
        ):
            min_size = (
                edata["type_def"].min_size
                or MoleculeType.types[Constants.SHARED_SYNTHESIZED_FRAGMENT].min_size
            )
            if edata["span"] > min_size * 2:
                new_edges += self.synthesis_partitioner(n1, n2, edata, border=min_size)

        syn_edges = [e for e in new_edges if e[2]["type_def"].design]
        AssemblyGraphBuilder.batch_add_edge_costs(
            self.graph, syn_edges, self.span_cost, None
        )

        for n1, n2, edata in syn_edges:
            self._complexity_to_efficiency(edata)

    def update(self):
        bad_edges = []
        for n1, n2, edata in self.graph.edges(data=True):
            self.update_long_pcr_products(n1, n2, edata)
            bad_edges += self.update_complexity(n1, n2, edata)

        self.logger.info(
            "Found {} highly complex synthesis segments".format(len(bad_edges))
        )

        # self.partition(bad_edges)

    # TODO: add logging to graph post processor
    # TODO: partition gaps
    def __call__(self):
        self.logger.info("Post processing graph for {}".format(self.query.name))
        self.update()
