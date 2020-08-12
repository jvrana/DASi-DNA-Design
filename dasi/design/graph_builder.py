import bisect
import itertools
from functools import partial
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union
from uuid import uuid4

import networkx as nx
import numpy as np
from Bio.SeqRecord import SeqRecord
from more_itertools import partition
from pyblast.utils import is_circular

from dasi.config import Config
from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.log import logger
from dasi.models import AlignmentContainer
from dasi.models import AlignmentGroup
from dasi.models import AssemblyNode
from dasi.models import MoleculeType
from dasi.models import MultiPCRProductAlignmentGroup
from dasi.models import PCRProductAlignmentGroup
from dasi.utils import bisect_between
from dasi.utils import Region
from dasi.utils import sort_with_keys
from dasi.utils.sequence import count_misprimings_in_amplicon
from dasi.utils.sequence import DNAStats
from dasi.utils.sequence.sequence_partitioner import find_by_partitions_for_sequence

SequenceScoringConfig = Config.SequenceScoringConfig


# TODO: whenever the efficiency is adjusted, record this in the notes
def add_edge_note(edata, key, value):
    if not edata["notes"]:
        edata["notes"] = {}
    edata["notes"][key] = value


Edge = Tuple[AssemblyNode, AssemblyNode, Dict]


class AssemblyGraphBuilder:
    """Class that builds an AssemblyGraph from an alignment container."""

    COST_THRESHOLD = Config.ASSEMBLY_COST_THRESHOLD

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

    def create_edge(
        self,
        n1: AssemblyNode,
        n2: AssemblyNode,
        cost: Union[float, None],
        material: Union[float, None],
        time: Union[float, None],
        efficiency: Union[float, None],
        span: int,
        type_def: MoleculeType,
        internal_or_external: str,
        condition: Tuple[bool, bool],
        group: Union[AlignmentGroup, PCRProductAlignmentGroup],
        notes: Dict = None,
        add_to_graph: bool = False,
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

        assert condition == type_def.design
        assert internal_or_external == type_def.int_or_ext
        edge = (
            n1,
            n2,
            dict(
                cost=cost,
                material=material,
                time=time,
                efficiency=efficiency,
                span=span,
                type_def=type_def,
                group=group,
                notes=notes,
            ),
        )
        if add_to_graph:
            self.G.add_edge(edge[0], edge[1], **edge[2])
        return edge

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
                            type_def=MoleculeType.types[group.type],
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
        """Add 'internal' edges from a list of AlignmentGroups. This means.

        :param groups:
        :return:
        """
        for g in groups:
            for a, b, ab_data in self.iter_internal_edge_data(g):
                for a_overhang, b_overhang in itertools.product(
                    [True, False], repeat=2
                ):
                    a_node = AssemblyNode(a[0], a[1], a[2], a_overhang)
                    b_node = AssemblyNode(b[0], b[1], b[2], b_overhang)
                    self.create_edge(a_node, b_node, add_to_graph=True, **ab_data)

    @staticmethod
    def make_overlap_iterator(
        a: Iterable[AssemblyNode], b: Iterable[AssemblyNode]
    ) -> Generator[Tuple[AssemblyNode, AssemblyNode], None, None]:
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

    @staticmethod
    def make_gap_iterator(
        a: Iterable[AssemblyNode], b: Iterable[AssemblyNode]
    ) -> Generator[Tuple[AssemblyNode, AssemblyNode], None, None]:
        """Find all nodes that satisfy the below condition: b.

        |--------|              |-------|              a
        """
        a, akeys = sort_with_keys(a, lambda x: x.index)
        b, bkeys = sort_with_keys(b, lambda x: x.index)

        for _a in a:
            i = bisect.bisect_right(bkeys, _a.index)
            for _b in b[:i]:
                yield _b, _a

    @staticmethod
    def make_origin_iterator(
        a: Iterable[AssemblyNode], b: Iterable[AssemblyNode], length: int
    ) -> Generator[Tuple[AssemblyNode, AssemblyNode], None, None]:
        a, akeys = sort_with_keys(a, lambda x: x.index)
        b, bkeys = sort_with_keys(b, lambda x: x.index)

        i = bisect.bisect_right(akeys, length)
        j = bisect.bisect_left(bkeys, length)
        for _a in a[:i]:
            for _b in b[j:]:
                yield _b, _a

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

        overlap_iter = self.make_overlap_iterator(a_nodes_overhang, b_nodes_overhang)
        gap_iter = self.make_gap_iterator(a_nodes_gap, b_nodes_gap)
        gap_origin_iter = self.make_origin_iterator(a_nodes_gap, b_nodes_gap, length)
        overlap_origin_iter = self.make_origin_iterator(
            a_nodes_overhang, b_nodes_overhang, length
        )

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
        self,
        bnode,
        anode,
        query_region,
        group_keys,
        groups,
        origin=False,
        validate_groups_present: bool = True,
        add_to_graph: bool = True,
    ):
        q = query_region.new(anode.index, bnode.index)
        if len(q) == 0:
            return
        filtered_groups = []
        if validate_groups_present:
            if not origin:
                i, j = bisect_between(group_keys, q.a, q.b)
                filtered_groups = groups[i:j]
            else:
                i = bisect.bisect_left(group_keys, q.a)
                j = bisect.bisect_right(group_keys, q.b)
                filtered_groups = groups[i:] + groups[:j]
        if not validate_groups_present or filtered_groups:
            if not origin:
                span = anode.index - bnode.index
            else:
                span = -len(q)
            if span <= 0:
                condition = (bnode.expandable, anode.expandable)
                return self.create_edge(
                    bnode,
                    anode,
                    cost=None,
                    material=None,
                    time=None,
                    efficiency=None,
                    internal_or_external="external",
                    type_def=MoleculeType.types[Constants.OVERLAP](condition),
                    condition=condition,
                    span=span,
                    group=None,
                    add_to_graph=add_to_graph,
                )

    def add_gap_edge(
        self, bnode, anode, query_region, origin=False, add_to_graph: bool = True
    ):
        """Add a 'gap' edge to the graph.

        :param bnode: ending AssemblyNode from the previous edge
        :param anode: starting AssemblyNode from the next edge
        :param query_region: a query region to copy (indices a,b does not matter)
        :param origin: whether this gap spans the origin or not
        :return:
        """
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
            return self.create_edge(
                bnode,
                anode,
                cost=None,
                material=None,
                time=None,
                efficiency=None,
                internal_or_external="external",
                type_def=MoleculeType.types[Constants.GAP](condition),
                span=span,
                condition=condition,
                group=None,
                add_to_graph=add_to_graph,
            )

    # TODO: refactor this. One method to update costs of edges. One method to
    #       add edges. One method to remove high cost edges.
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
                    if _d["material"] is None:
                        raise Exception
                    e.update(_d)
        if cost_threshold is not None:
            graph.remove_edges_from(edges_to_remove)

    def _get_cost(self, bp, ext):
        data = self.span_cost(bp, ext).data
        return {k: v[0] for k, v in data.items()}

    def update_costs(self, edges: Optional[List[Edge]] = None):
        if edges is None:
            edges = []
            for n1, n2, edata in self.G.edges(data=True):
                if edata["cost"] is None:
                    edges.append((n1, n2, edata))
                elif edata["type_def"].name == Constants.SHARED_SYNTHESIZED_FRAGMENT:
                    edges.append((n1, n2, edata))
        self.batch_add_edge_costs(self.G, edges, self.span_cost, self.COST_THRESHOLD)

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
        self.update_costs()

        # TODO: freeze?
        # nx.freeze(self.G)
        return self.G


# TODO: evaluate primer designs, scoring PCR products d(how long does this take?)
# TODO: refactor this class, expose config options


class AssemblyGraphPostProcessor:
    """Pre-processing for assembly graphs. Evaluates:

    1. synthesis complexity and weights corresponding edge
    2. pcr product efficiency
    3. (optional) optimal partitions for synthesis fragments
    """

    # TODO: add post processing config

    def __init__(
        self,
        graph: nx.DiGraph,
        query: SeqRecord,
        span_cost: SpanCost,
        seqdb: Dict[str, SeqRecord],
        container: AlignmentContainer,
        stats_repeat_window: Optional[int] = None,
        stats_window: Optional[int] = None,
        stats_hairpin_window: Optional[int] = None,
        edge_threshold: Optional[float] = None,
        stages: Optional[Tuple[str]] = None,
    ):
        if stats_repeat_window is None:
            stats_repeat_window = SequenceScoringConfig.stats_repeat_window
        if stats_window is None:
            stats_window = SequenceScoringConfig.stats_window
        if stats_hairpin_window is None:
            stats_hairpin_window = SequenceScoringConfig.stats_hairpin_window
        if stages is None:
            stages = SequenceScoringConfig.post_process_stages
        if edge_threshold is None:
            edge_threshold = SequenceScoringConfig.edge_threshold
        self.graph = graph
        self.graph_builder = AssemblyGraphBuilder(container, span_cost=span_cost)
        self.graph_builder.G = graph
        self.query = query
        self.seqdb = seqdb
        query_seq = str(query.seq)
        if is_circular(query):
            query_seq = query_seq + query_seq
        self.stats = DNAStats(
            query_seq,
            repeat_window=stats_repeat_window,
            stats_window=stats_window,
            hairpin_window=stats_hairpin_window,
        )
        self.stats_single = DNAStats(
            str(query.seq),
            repeat_window=stats_repeat_window,
            stats_window=stats_window,
            hairpin_window=stats_hairpin_window,
        )
        self.logged_msgs = []
        # TODO: make a more sophisticated complexity function?
        # TODO: expose this to input parameters
        self.COMPLEXITY_THRESHOLD = SequenceScoringConfig.complexity_threshold
        self.logger = logger(self)
        self.span_cost = span_cost
        self.stages = stages
        self.edge_threshold = edge_threshold

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

    @staticmethod
    def _adj_eff(edata, e):
        edata["efficiency"] = e
        if e == 0:
            edata["cost"] = np.inf
        if edata["material"] is None:
            edata["cost"] = np.inf
        else:
            edata["cost"] = edata["material"] / edata["efficiency"]

    def _complexity_to_efficiency(self, edata: dict) -> bool:
        if not edata["complexity"] >= 0.0:
            raise ValueError("Complexity is not defined. {}".format(edata))
        ratio = edata["complexity"] / self.COMPLEXITY_THRESHOLD
        if ratio >= 1.0:
            e = SequenceScoringConfig.not_synthesizable_efficiency / ratio
            self._adj_eff(edata, e)
            return True
        return False

    @staticmethod
    def _is_pcr_product(edata):
        return edata["type_def"].name in [
            Constants.PCR_PRODUCT,
            Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
            Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
            Constants.PCR_PRODUCT_WITH_PRIMERS,
        ]

    def score_long_pcr_products(self, n1, n2, edata):
        if self._is_pcr_product(edata):
            span = edata["span"]
            for a, b, c in SequenceScoringConfig.pcr_length_range_efficiency_multiplier:
                if a <= span < b:
                    self._adj_eff(edata, edata["efficiency"] * c)
                    add_edge_note(edata, "long_pcr_product", True)
                    break

    def _score_misprimings_from_alignment(self, alignment):
        subject_key = alignment.subject_key
        subject = self.seqdb[subject_key]
        i = alignment.subject_region.a
        j = alignment.subject_region.b
        subject_seq = str(subject.seq)

        return count_misprimings_in_amplicon(
            subject_seq,
            i,
            j,
            min_primer_anneal=SequenceScoringConfig.mispriming_min_anneal,
            max_primer_anneal=SequenceScoringConfig.mispriming_max_anneal,
            cyclic=alignment.subject_region.cyclic,
        )

    # TODO: select the best template...
    # TODO: speed up this process
    def score_primer_misprimings(self, n1, n2, edata):
        if self._is_pcr_product(edata):

            group = edata["group"]

            misprime_list = []

            if isinstance(group, MultiPCRProductAlignmentGroup):
                prioritize_function = group.prioritize_groupings
                alignments = group.iter_templates()

            elif isinstance(group, AlignmentGroup):
                prioritize_function = group.prioritize_alignments
                alignments = group.alignments
            else:
                raise TypeError(
                    "Group '{}' not supported by this function.".format(group.__class__)
                )

            arr = []
            for index, alignment in enumerate(alignments):
                if not (
                    "PCR" in alignment.type or alignment.type == Constants.FRAGMENT
                ):
                    raise ValueError(
                        "{} is not a valid type for a template".format(alignment.type)
                    )
                mispriming = self._score_misprimings_from_alignment(alignment)
                arr.append((mispriming, index, alignment))
                if mispriming == 0:
                    break
            arr.sort(key=lambda x: x[0])
            indices = [a[1] for a in arr]
            prioritize_function(indices)
            score = arr[0][0]

            self._adj_eff(
                edata,
                edata["efficiency"] * SequenceScoringConfig.mispriming_penalty ** score,
            )
            add_edge_note(edata, "num_misprimings", score)
            add_edge_note(edata, "n_templates_eval", len(misprime_list))

    def _filter_partition_edges(self, edges: List[Edge]) -> List[Edge]:
        edges_to_partition = []
        for n1, n2, edata in edges:
            min_size = (
                edata["type_def"].min_size
                or MoleculeType.types[Constants.SHARED_SYNTHESIZED_FRAGMENT].min_size
            )
            if (
                edata["span"] > min_size * 2
                and edata["complexity"]
                >= Config.SequenceScoringConfig.complexity_threshold
            ):
                edges_to_partition.append((n1, n2, edata))
        return edges_to_partition

    def partition(self, edges: List[Edge]):
        tracker = self.logger.track(
            "INFO", desc="Partitioning sequences", total=3
        ).enter()
        tracker.update(0, "{} highly complex sequences".format(len(edges)))

        edges_to_partition = self._filter_partition_edges(edges)

        cyclic = is_circular(self.query)
        partitions = find_by_partitions_for_sequence(
            self.stats_single,
            cyclic=cyclic,
            threshold=Config.SequenceScoringConfig.complexity_threshold,
            step_size=Config.SequenceScoringConfig.partition_step_size,
            delta=Config.SequenceScoringConfig.partition_overlap,
        )
        tracker.update(1, "Partition: locations: {}".format(partitions))
        add_gap_edge = partial(self.graph_builder.add_gap_edge, add_to_graph=False)
        add_overlap_edge = partial(
            self.graph_builder.add_overlap_edge,
            add_to_graph=True,
            validate_groups_present=False,
            groups=None,
            group_keys=None,
        )

        new_edges = []
        for n1, n2, edata in edges_to_partition:
            r = Region(n1.index, n2.index, len(self.query.seq), cyclic=cyclic)

            if n1.type == "B" and n2.type == "A":

                for p in partitions:
                    if p in r:
                        # TODO: overlap? find optimal partition for overlap?
                        i4 = p
                        i3 = p + Config.SequenceScoringConfig.partition_overlap

                        n3 = AssemblyNode(i3, False, "BC", overhang=True)
                        n4 = AssemblyNode(i4, False, "CA", overhang=True)
                        e1 = add_gap_edge(n1, n3, r, origin=False)
                        e2 = add_gap_edge(n1, n3, r, origin=True)
                        if e1 is None and e2 is None:
                            continue
                        e3 = add_overlap_edge(n3, n4, r, origin=False)
                        e4 = add_overlap_edge(n3, n4, r, origin=True)
                        if e3 is None and e4 is None:
                            continue
                        e5 = add_gap_edge(n4, n2, r, origin=False)
                        e6 = add_gap_edge(n4, n2, r, origin=True)
                        if e5 is None and e6 is None:
                            continue
                        new_edges += [e1, e2, e3, e4, e5, e6]
        for e in new_edges:
            if e is not None:
                self.graph_builder.G.add_edge(e[0], e[1], **e[2])
        edges = []
        for n1, n2, edata in self.graph_builder.G.edges(data=True):
            if edata["cost"] is None:
                edges.append((n1, n2, edata))
        tracker.update(2, "Partition: Added {} new edges".format(len(edges)))
        self.graph_builder.update_costs(edges)
        self.score_complexity_edges(list(self.graph_builder.G.edges(data=True)))
        tracker.exit()

    def remove_inefficient_edges(self):
        to_remove = []
        for n1, n2, edata in self.graph_builder.G.edges(data=True):
            if edata["efficiency"] < self.edge_threshold:
                to_remove.append((n1, n2))
        self.graph_builder.G.remove_edges_from(to_remove)
        self.logger.info("Removed {} inefficient edges".format(len(to_remove)))

    def _score_syn_dna(
        self, n1: AssemblyNode, n2: AssemblyNode, edata: dict
    ) -> List[Edge]:
        if edata["type_def"].synthesize:
            span = edata["span"]
            if span > 0:
                # TODO: cyclic may not always be true
                region = self._edge_to_region(n1, n2)
                a = region.a
                c = region.c
                if c < a:
                    c += region.context_length
                assert c <= len(self.stats)
                complexity = self.stats.cost(a, c)
                if not complexity >= 0.0:
                    pass
                edata["complexity"] = complexity
                if self._complexity_to_efficiency(edata):
                    add_edge_note(edata, "highly_complex", True)
                    return True
        return False

    def score_complexity_edges(self, edges: List[Edge] = None) -> List[Edge]:
        """Score synthetic edges."""
        if edges is None:
            edges = self.graph.edges(data=True)
        complex_edges = []
        for n1, n2, edata in self.logger.tqdm(
            edges, "INFO", desc="Scoring synthetic DNA"
        ):
            if self._score_syn_dna(n1, n2, edata):
                complex_edges.append((n1, n2, edata))
        return complex_edges

    def update(self):
        self.logger.info("Post Processor: {}".format(self.stages))
        bad_edges = []

        self.logger.info("Scoring long PCR products")
        edges = list(self.graph.edges(data=True))

        if SequenceScoringConfig.SCORE_LONG in self.stages:
            for n1, n2, edata in self.logger.tqdm(
                edges, "INFO", desc="Scoring long PCR products"
            ):
                self.score_long_pcr_products(n1, n2, edata)
        if SequenceScoringConfig.SCORE_MISPRIMINGS in self.stages:
            for n1, n2, edata in self.logger.tqdm(
                edges, "INFO", desc="Scoring primer misprimings"
            ):
                self.score_primer_misprimings(n1, n2, edata)
        if SequenceScoringConfig.SCORE_COMPLEXITY in self.stages:
            bad_edges += self.score_complexity_edges(edges)
        self.logger.info(
            "Found {} highly complex synthesis segments".format(len(bad_edges))
        )

        if SequenceScoringConfig.PARTITION in self.stages:
            self.partition(bad_edges)

        # TODO: reimplement `remove_inefficient_edges`
        # self.remove_inefficient_edges()

    # TODO: add logging to graph post processor
    # TODO: partition gaps
    def __call__(self):
        self.logger.info("Post processing graph for {}".format(self.query.name))
        self.update()
