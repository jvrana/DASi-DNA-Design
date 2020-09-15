from copy import deepcopy
from typing import Dict
from typing import List
from typing import Tuple

import networkx as nx

from dasi.constants import Constants
from dasi.design import DesignResult
from dasi.design.design import DesignABC
from dasi.design.graph_builder import AssemblyGraphPostProcessor
from dasi.models import Alignment
from dasi.models import AlignmentContainer
from dasi.models import AlignmentGroup
from dasi.utils import group_by
from dasi.utils import log_metadata
from dasi.utils import sort_with_keys


def overlapping_groups(group_list_a, group_list_b):
    """Get all groups in group_list_b that right-hand overlap with
    group_list_a."""
    group_sort, group_keys = sort_with_keys(
        group_list_b, key=lambda x: x.query_region.a
    )
    tuples = []
    for group_a in group_list_a:
        overlapping = AlignmentContainer.filter_alignments_by_span(
            group_sort,
            group_a.query_region,
            key=lambda p: p.query_region.a,
            end_inclusive=False,
        )

        if group_a in overlapping:
            overlapping.remove(group_a)
        tuples.append((group_a, overlapping))
    return tuples


def add_clusters(design):
    # list of all alignment groups
    all_groups = []
    for container in design.container_factory.containers().values():
        all_groups += container.get_groups_by_types(Constants.SHARED_FRAGMENT)

    # alignment_groups grouped by query_key
    grouped_by_qk = {}
    for g in all_groups:
        grouped_by_qk.setdefault(g.query_key, list())
        grouped_by_qk[g.query_key].append(g)

    # overlapping_by_qk
    overlapping = []
    for qk, groups in list(grouped_by_qk.items())[:]:
        overlapping += overlapping_groups(groups, groups)

    new_alignments = []
    for group_a, group_list in overlapping:
        new = container.expand_overlaps(
            group_list + [group_a], include_left=False, atype=Constants.SHARED_FRAGMENT
        )
        new_alignments += new
    new_alignments = list({a.eq_hash(): a for a in new_alignments}.values())
    design.container_factory.add_alignments(new_alignments)


def to_undirected(graph):
    """.to_undirected is implemented in networkx out of the box, however, it
    suffers from occational infinite recursion errors during the deepcopy phase
    of the method (unknown as to why)."""
    undirected = nx.Graph()
    copied = deepcopy(graph)
    for n in copied.nodes:
        ndata = copied.nodes[n]
        undirected.add_node(n, **ndata)
    for n1, n2 in copied.edges:
        edata = copied.edges[n1, n2]
        undirected.add_edge(n1, n2, **edata)
    return undirected


def get_subgraphs(graph):
    """Get independent subgraphs."""
    node_list = list(graph.nodes)
    subgraphs = []
    while len(node_list) > 0:
        node = node_list[-1]
        subgraph = nx.bfs_tree(to_undirected(graph), node)
        for n in subgraph.nodes:
            node_list.remove(n)
        subgraphs.append(graph.subgraph(subgraph.nodes))
    return subgraphs


def has_repeats(g):
    """Check if the interaction graph has a repeated DNA sequence."""
    grouped_by_key = {}
    for n in g.nodes:
        grouped_by_key.setdefault(n[0], list())
        grouped_by_key[n[0]].append((n[1], n[2]))
    for k, v in grouped_by_key.items():
        if len(v) > 1:
            print(grouped_by_key)
            return True
    return False


def cluster_graph(design):
    interaction_graph = nx.Graph()
    all_groups = []
    for container in design.container_factory.containers().values():
        all_groups += container.get_groups_by_types(Constants.SHARED_FRAGMENT)
    for g in all_groups:
        for a in g.alignments:
            n1 = (a.query_key, a.query_region.a, a.query_region.b, a.query_region.c)
            n2 = (
                a.subject_key,
                a.subject_region.a,
                a.subject_region.b,
                a.subject_region.c,
            )
            edata = interaction_graph.get_edge_data(n1, n2)
            if edata:
                edata.setdefault("alignments", list())
                edata["alignments"].append(a)
            else:
                interaction_graph.add_edge(n1, n2, alignments=[a])
    graphs = get_subgraphs(interaction_graph)
    graphs = [g for g in graphs if not has_repeats(g)]
    graphs.sort(reverse=True, key=lambda x: x.number_of_nodes())
    return graphs


class LibraryDesign(DesignABC):
    """Design class for producing assemblies for libraries."""

    # TODO: move to config
    FAVOR_SHARED_SEQUENCES = 2
    DEFAULT_N_ASSEMBLIES = DesignABC.DEFAULT_N_ASSEMBLIES
    ALGORITHM = Constants.ALGORITHM_LIBRARY

    def __init__(self, span_cost=None, seqdb=None):
        """Something."""
        super().__init__(span_cost=span_cost, seqdb=seqdb)
        self.shared_alignments = []
        self._edges = []

    def _expand_from_synthesized(self):
        """Expand PCR products that are next to SYNTHESIZED_FRAGMENTS.

        :return:
        """
        for query_key, container in self.container_factory.containers().items():

            def is_not_shared(group_a: AlignmentGroup, group_b: AlignmentGroup):
                if group_a.type == Constants.SHARED_SYNTHESIZED_FRAGMENT:
                    return False
                if group_b.type == Constants.SHARED_SYNTHESIZED_FRAGMENT:
                    return True
                return False

            new_alignments = container.expand_overlaps(
                container.get_groups_by_types(
                    [
                        Constants.FRAGMENT,
                        Constants.PCR_PRODUCT,
                        Constants.SHARED_SYNTHESIZED_FRAGMENT,
                    ]
                ),
                Constants.PCR_PRODUCT,
                pass_condition=is_not_shared,
            )
            container.add_alignments(new_alignments, lim_size=True)

            self.logger.info(
                "{}: Expanded {} using {} and found {} new alignments.".format(
                    query_key,
                    Constants.PCR_PRODUCT,
                    Constants.SHARED_SYNTHESIZED_FRAGMENT,
                    len(new_alignments),
                )
            )

            # # grab the pcr products and expand primer pairs (again)
            # templates = container.get_groups_by_types(Constants.PCR_PRODUCT)
            # new_primer_pairs = container.expand_primer_pairs(templates)
            # self.logger.info(
            #     "{}: Expanded {} {} using {}".format(
            #         query_key,
            #         len(new_primer_pairs),
            #         Constants.PCR_PRODUCT,
            #         Constants.SYNTHESIZED_FRAGMENT,
            #     )
            # )

    def _expand_synthesized_fragments(self):
        """
        1. copy groups from SHARED_FRAGMENTS to SYNTHESIZED_FRAGMENT
        2. expand overlaps for SYNTHESIZED_FRAGMENTS
        """
        for query_key, container in self.container_factory.containers().items():
            # expand the share fragments using their own endpoints
            original_shared_fragments = container.get_groups_by_types(
                Constants.SHARED_FRAGMENT
            )
            copied_alignments = container.copy_groups(
                original_shared_fragments, Constants.SHARED_SYNTHESIZED_FRAGMENT
            )
            container.add_alignments(copied_alignments, lim_size=True)

            new_alignments = container.expand_overlaps(
                container.get_groups_by_types(Constants.SHARED_SYNTHESIZED_FRAGMENT),
                atype=Constants.SHARED_SYNTHESIZED_FRAGMENT,
            )

            self.logger.info(
                "{}: Expanded {} new {} alignments.".format(
                    query_key,
                    len(new_alignments),
                    Constants.SHARED_SYNTHESIZED_FRAGMENT,
                )
            )
            container.add_alignments(new_alignments, lim_size=True)

    def _check_shared_repeats(self):
        repeats = []
        for query_key, container in self.container_factory.containers().items():
            for align in container.get_alignments_by_types(Constants.SHARED_FRAGMENT):
                qk = align.query_key
                sk = align.subject_key
                if qk == sk:
                    repeats.append(align)
        assert not repeats

    def _share_query_blast(self):
        """Find and use shared fragments across queries.

        :return:
        """

        # step 1: get query-on-query alignments
        self.logger.info("=== Expanding shared library fragments ===")
        blast = self.blast_factory(self.QUERIES, self.QUERIES)
        blast.blastn()

        results = blast.get_perfect()

        # step 2: eliminate self binding results
        results = [
            entry
            for entry in results
            if entry["query"]["origin_key"] != entry["subject"]["origin_key"]
        ]

        self.logger.info(
            "Found {} shared alignments between the queries".format(len(results))
        )

        # step 3: load results to the container
        self.shared_alignments = results
        self.container_factory.seqdb.update(blast.seq_db.records)
        self.container_factory.load_blast_json(results, Constants.SHARED_FRAGMENT)

    def _precompile_library_expansion(self):
        self._share_query_blast()
        self._expand_synthesized_fragments()
        self._expand_from_synthesized()
        self._check_shared_repeats()

    def precompile(self):
        self.uncompile()
        tracker = self.logger.track(
            "INFO", desc="Precompiling library", total=5
        ).enter()

        self.graphs = {}

        tracker.update(0, "Running blast")
        self._blast()

        tracker.update(1, "Running shared fragment blast")
        self._share_query_blast()

        # tracker.update(2, "Expanding shared fragments")
        # self._expand_synthesized_fragments()
        self._expand_from_synthesized()

        tracker.update(3, "Finding shared clusters")
        self.update_library_metadata()
        tracker.exit()

    def _adjust_shared_synthetic_fragment(self):
        for qk, graph in self.graphs.items():
            for n1, n2, edata in graph.edges(data=True):
                if edata["type_def"].name == Constants.SHARED_SYNTHESIZED_FRAGMENT:
                    group = edata["group"]

                    # TODO: better way to make notes on graph edge
                    if "notes" not in edata or edata["notes"] is None:
                        edata["notes"] = {}

                    # if "n_clusters" not in group.meta:
                    #     pass
                    edata["notes"]["n_clusters"] = group.meta["n_clusters"]
                    # TODO: adjust n_clusters
                    edata["material"] = (
                        edata["material"]
                        / (group.meta["n_clusters"])
                        / self.FAVOR_SHARED_SEQUENCES
                    )
                    edata["cost"] = edata["material"] / edata["efficiency"]

    def postcompile(self, post_compile_kwargs: dict = None):
        self.logger.info("Post-processing graphs")
        self.post_process_graphs(post_compile_kwargs)
        self._adjust_shared_synthetic_fragment()

    def update_library_metadata(self):
        add_clusters(self)
        graphs = cluster_graph(self)

        for c in self.container_list:
            c.share_group_tag = {}

        # update the meta data
        for g in graphs:
            n_clusters = g.number_of_nodes()
            for n1, n2, edata in g.edges(data=True):
                alignments = edata["alignments"]
                for qk, container_alignments in group_by(
                    alignments, key=lambda x: x.query_key
                ).items():
                    container = self.containers[qk]
                    for a in container_alignments:
                        group_key = (
                            a.query_region.a,
                            a.query_region.b,
                            Constants.SHARED_SYNTHESIZED_FRAGMENT,
                        )
                        container.group_tags.add(
                            Constants.SHARE_GROUP_TAG,
                            group_key,
                            {
                                "alignments": container_alignments,
                                "cross_container_alignments": alignments,
                                "n_clusters": n_clusters,
                            },
                        )

        # for container in design.container_list():
        #     groups = container.get_groups_by_types(Constants.SHARED_FRAGMENT)
        #     copied = container.copy_groups(groups, atype=Constants.SYNTHESIZED_FRAGMENT)
        #     container.add_alignments(copied)

        # def is_not_shared(group_a: AlignmentGroup, group_b: AlignmentGroup):
        #     if group_a.type == Constants.SHARED_FRAGMENT:
        #         return False
        #     if group_b.type == Constants.SYNTHESIZED_FRAGMENT:
        #         return True
        #     return False

        # new_alignments = container.expand_overlaps(
        #     container.get_groups_by_types(
        #         [
        #             Constants.FRAGMENT,
        #             Constants.PCR_PRODUCT,
        #             Constants.SYNTHESIZED_FRAGMENT,
        #         ]
        #     ),
        #     Constants.PCR_PRODUCT,
        #     pass_condition=is_not_shared,
        # )
        # container.add_alignments(new_alignments, lim_size=True)

    # def post_process_library_graphs(self):
    #     for qk, graph in self.graphs.items():
    #         query = self.seqdb[qk]
    #         processor = AssemblyGraphPostProcessor(graph, query)
    #         processor()

    # @log_metadata(
    #     "optimize", additional_metadata={"algorithm": Constants.ALGORITHM_LIBRARY}
    # )
    # def optimize(
    #     self, n_paths: int = DEFAULT_N_ASSEMBLIES, n_jobs: int = DEFAULT_N_JOBS
    # ) -> Dict[str, DesignResult]:
    #     """Optimize the assembly graph for library assembly."""
    #     return super().optimize(n_paths=n_paths, n_jobs=n_jobs)
