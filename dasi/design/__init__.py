"""Primer and synthesis design

.. module:: design

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    assembly
"""
from __future__ import annotations
from dasi.alignments import (
    Alignment,
    AlignmentContainerFactory,
    AlignmentContainer,
    AlignmentGroup,
    ComplexAlignmentGroup,
)
from dasi.constants import Constants
from .graph_builder import AssemblyGraphBuilder
from dasi.utils import (
    perfect_subject,
    multipoint_shortest_path,
    sort_with_keys,
    sort_cycle,
)
from dasi.exceptions import DasiDesignException
import networkx as nx
from pyblast import BioBlastFactory
from dasi.log import logger
from typing import List, Tuple, Dict
from Bio.SeqRecord import SeqRecord
import numpy as np
from more_itertools import pairwise
from pyblast.utils import Span, is_circular
import pandas as pd
import bisect
from dasi.design.graph_builder import AssemblyNode
from collections.abc import Iterable
from itertools import zip_longest
from copy import deepcopy


BLAST_PENALTY_CONFIG = {"gapopen": 3, "gapextend": 3, "reward": 1, "penalty": -5}


class DesignResult(Iterable):
    def __init__(self, container, graph, query_key):
        self.container = container
        self.graph = graph
        self.query_key = query_key
        self.query = self.container.seqdb[query_key]
        self._assemblies = []
        self._keys = []

    @property
    def assemblies(self):
        return tuple(self._assemblies)

    def _new(self, path: List[AssemblyNode]):
        return Assembly(path, self.container, self.graph, self.query_key, self.query)

    def add_assembly(self, path: List[AssemblyNode]):
        assembly = self._new(path)
        cost = assembly.cost()
        n_nodes = len(assembly._nodes)
        k = (cost, n_nodes)
        i = bisect.bisect_left(self._keys, k)
        self._assemblies.insert(i, assembly)
        self._keys.insert(i, k)

    def add_assemblies(self, paths: List[List[AssemblyNode]]):
        for path in paths:
            self.add_assembly(path)

    def __iter__(self):
        for assembly in self.assemblies:
            yield assembly

    def __getitem__(self, item):
        return list(self)[item]


class Assembly(Iterable):
    """
    Should take in a path, graph, container, seqdb to produce relevant information
    """

    def __init__(
        self,
        nodes: List[AssemblyNode],
        container: AlignmentContainer,
        full_assembly_graph: nx.DiGraph,
        query_key: str,
        query: SeqRecord,
    ):
        self.container = container
        self.query_key = query_key
        self.query = query
        self._nodes = tuple(nodes)
        self._full_graph = full_assembly_graph
        self.graph = self._subgraph(self._full_graph, nodes)
        nx.freeze(self.graph)

    def _subgraph(self, graph: nx.DiGraph, nodes: List[AssemblyNode]):
        def _resolve(node: AssemblyNode, query_region) -> Tuple[int, dict]:
            new_node = AssemblyNode(query_region.t(node.index), *list(node)[1:])
            return new_node, {}

        SG = nx.OrderedDiGraph()
        nodes = [AssemblyNode(*n) for n in nodes]
        example_query_region = self.container.alignments[0].query_region

        resolved_nodes = [_resolve(node, example_query_region) for node in nodes]
        if self.cyclic:
            resolved_nodes = sort_cycle(
                resolved_nodes, key=lambda n: (n[0].type, n[0].index, n)
            )
        SG.add_nodes_from(resolved_nodes)

        if self.cyclic:
            pair_iter = list(pairwise(nodes + nodes[:1]))
        else:
            pair_iter = list(pairwise(nodes))

        for n1, n2 in pair_iter:
            edata = deepcopy(graph.get_edge_data(n1, n2))
            if edata is None:
                edata = {
                    "weight": np.inf,
                    "type": "missing",
                    "span": np.inf,
                    "name": "missing",
                }

            # TODO: fix query_region (overlaps are backwards)
            query_region = self.container.alignments[0].query_region.new(
                n1.index, n2.index, allow_wrap=True
            )
            groups = self.container.find_groups_by_pos(query_region.a, query_region.b)
            edata["groups"] = groups
            edata["query_region"] = query_region
            SG.add_edge(
                _resolve(n1, query_region)[0], _resolve(n2, query_region)[0], **edata
            )
        return SG

    @property
    def cyclic(self):
        return is_circular(self.query)

    def cost(self):
        total = 0
        for _, _, edata in self.edges():
            total += edata["weight"]
        return total

    def edges(self, data=True) -> Iterable[Tuple[AssemblyNode, AssemblyNode, Dict]]:
        return self.graph.edges(data=data)

    def nodes(self, data=True) -> Iterable[Tuple[AssemblyNode, Dict]]:
        return self.graph.nodes(data=data)

    def edit_distance(
        self, other: Assembly, explain=False
    ) -> Tuple[int, List[Tuple[int, str]]]:
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
                if isinstance(group, ComplexAlignmentGroup):
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
                    "weight": edata["weight"],
                    "span": edata["span"],
                    "type": edata["type"],
                    "name": edata["name"],
                }
            )

        df = pd.DataFrame(rows)
        return df

    def __eq__(self, other: Assembly) -> bool:
        return self.edit_distance(other) == 0

    def __iter__(self):
        for n in self.nodes(data=False):
            yield n


class Design(object):
    """
    Design class that returns optimal assemblies from a set of materials.
    """

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"
    FRAGMENTS = "fragments"

    def __init__(self, span_cost=None, seqdb=None):
        self.blast_factory = BioBlastFactory()
        self.logger = logger(self)

        # graph by query_key
        if seqdb is None:
            seqdb = {}
        self._seqdb = seqdb
        self.span_cost = span_cost
        self.graphs = {}
        self.results = {}
        self.container_factory = AlignmentContainerFactory(self.seqdb)

    @property
    def seqdb(self):
        return self._seqdb

    def add_materials(
        self,
        primers: List[SeqRecord],
        templates: List[SeqRecord],
        queries: List[SeqRecord],
        fragments=None,
    ):
        if fragments is None:
            fragments = []
        self.add_primers(primers)
        fragments = self.filter_linear_records(fragments)
        self.add_templates(templates + fragments)
        self.add_queries(queries)
        self.add_fragments(fragments)

        self.template_results = []
        self.fragment_results = []
        self.primer_results = []

    def add_primers(self, primers: List[SeqRecord]):
        self.logger.info("Adding primers")
        self.blast_factory.add_records(primers, self.PRIMERS)

    def add_templates(self, templates: List[SeqRecord]):
        self.logger.info("Adding templates")
        self.blast_factory.add_records(templates, self.TEMPLATES)

    def add_queries(self, queries: List[SeqRecord]):
        self.logger.info("Adding queries")
        self.blast_factory.add_records(queries, self.QUERIES)

    def add_fragments(self, fragments: List[SeqRecord]):
        self.logger.info("Adding fragments")
        self.blast_factory.add_records(fragments, self.FRAGMENTS)

    @classmethod
    def filter_linear_records(cls, records):
        """Return only linear records"""
        return [r for r in records if not is_circular(r)]

    @classmethod
    def filter_perfect_subject(cls, results):
        """return only results whose subject is 100% aligned to query"""
        return [r for r in results if perfect_subject(r["subject"])]

    # TODO: do a single blast and sort results based on record keys
    def _blast(self):
        self.logger.info("Compiling assembly graph")

        # align templates
        blast = self.blast_factory(self.TEMPLATES, self.QUERIES)
        blast.update_config(BLAST_PENALTY_CONFIG)
        blast.quick_blastn()
        results = blast.get_perfect()
        self.template_results = results

        # align fragments
        if self.blast_factory.record_groups[self.FRAGMENTS]:
            fragment_blast = self.blast_factory(self.FRAGMENTS, self.QUERIES)
            fragment_blast.update_config(BLAST_PENALTY_CONFIG)
            fragment_blast.quick_blastn()
            fragment_results = blast.get_perfect()
            fragment_results = self.filter_perfect_subject(fragment_results)
        else:
            fragment_results = []
        self.fragment_results = fragment_results

        self.container_factory.seqdb.update(blast.seq_db.records)
        self.logger.info("Number of template matches: {}".format(len(results)))
        self.logger.info(
            "Number of perfect fragment matches: {}".format(len(fragment_results))
        )

        # align primers
        if self.blast_factory.record_groups[self.PRIMERS]:
            primer_blast = self.blast_factory(self.PRIMERS, self.QUERIES)
            primer_blast.update_config(BLAST_PENALTY_CONFIG)
            primer_blast.quick_blastn_short()
            primer_results = primer_blast.get_perfect()
            primer_results = self.filter_perfect_subject(primer_results)
            self.container_factory.seqdb.update(primer_blast.seq_db.records)
            self.logger.info(
                "Number of perfect primers: {}".format(len(primer_results))
            )
        else:
            primer_results = []
        self.primer_results = primer_results

        self.container_factory.load_blast_json(fragment_results, Constants.FRAGMENT)
        self.container_factory.load_blast_json(results, Constants.PCR_PRODUCT)
        self.container_factory.load_blast_json(primer_results, Constants.PRIMER)

    @property
    def containers(self):
        return self.container_factory.containers()

    def container_list(self):
        return list(self.container_factory.containers().values())

    def query_keys(self):
        return list(self.container_factory.containers())

    def assemble_graphs(self):
        for query_key, container in self.logger.tqdm(
            self.container_factory.containers().items(),
            "INFO",
            desc="compiling all containers",
        ):
            container.expand(expand_overlaps=True, expand_primers=True)

            # group by query_regions
            groups = container.groups()

            self.logger.info(
                "Number of types: {}".format(len(container.groups_by_type))
            )
            self.logger.info("Number of groups: {}".format(len(groups)))

            # build assembly graph
            graph_builder = AssemblyGraphBuilder(container, span_cost=self.span_cost)
            G = graph_builder.build_assembly_graph()

            self.logger.info("=== Assembly Graph ===")
            self.logger.info(nx.info(G))
            assert G.number_of_edges()
            self.graphs[query_key] = G

    def compile(self):
        """Compile materials to assembly graph"""
        self.results = {}
        self._blast()
        self.assemble_graphs()

    # def plot_matrix(self, matrix):
    # plot matrix
    # import pylab as plt
    # import seaborn as sns
    # import numpy as np
    #
    # plot_matrix = matrix.copy()
    # plot_matrix[plot_matrix == np.inf] = 10000
    # plot_matrix = np.nan_to_num(plot_matrix)
    #
    # fig = plt.figure(figsize=(24, 20))
    # ax = fig.gca()
    # step = 1
    # sns.heatmap(plot_matrix[::step, ::step], ax=ax)

    @staticmethod
    def _find_iter_alignment(a, b, alignments):
        for align in alignments:
            if a == align.query_region.a and b == align.query_region.b:
                yield align

    def _fragment(self, query_key, a, b, fragment_type, cost):
        def sub_record(record, span):
            ranges = span.ranges()
            sub = record[ranges[0][0] : ranges[0][1]]
            for r in ranges[1:]:
                sub += record[r[0] : r[1]]
            sub.annotations = record.annotations
            return sub

        alignments = self.container_factory.alignments[query_key]
        align = list(self._find_iter_alignment(a, b, alignments))[0]
        subject_key = align.subject_key
        subject_rec = self.container_factory.seqdb[subject_key]
        query_rec = self.container_factory.seqdb[query_key]

        subject_seq = sub_record(subject_rec, align.subject_region)

        fragment_info = {
            "query_id": query_key,
            "query_name": query_rec.name,
            "query_region": (align.query_region.a, align.query_region.b),
            "subject_id": subject_key,
            "subject_name": subject_rec.name,
            "subject_region": (align.subject_region.a, align.subject_region.b),
            "fragment_length": len(align.subject_region),
            "fragment_seq": subject_seq,
            "cost": cost,
            "type": fragment_type,
        }

    def path_to_edge_costs(self, path, graph):
        arr = []
        for n1, n2 in pairwise(path):
            edata = graph[n1][n2]
            arr.append((n1, n2, edata))
        return arr

    def path_to_df(self, paths_dict):

        fragments = []
        primers = []

        for qk, paths in paths_dict.items():
            paths = paths
            G = self.graphs[qk]
            container = self.container_factory.containers()[qk]

            record = self.container_factory.seqdb[qk]
            path = paths[0] + paths[0][:1]

            for n1, n2 in pairwise(path):
                edata = G[n1][n2]
                cost = edata["weight"]
                if n1[2] == "A" and n2[2] == "B":
                    A = n1[0]
                    B = n2[0]
                    group = container.find_groups_by_pos(A, B)[0]

                    if isinstance(group, AlignmentGroup):
                        align = group.alignments[0]
                        subject = align.subject_key
                        subject_rec = self.container_factory.seqdb[align.subject_key]
                        subject_rec_name = subject_rec.name
                        subject_seq = str(
                            subject_rec[
                                align.subject_region.a : align.subject_region.b
                            ].seq
                        )
                        subject_region = (
                            align.subject_region.a,
                            align.subject_region.b,
                        )
                    elif isinstance(group, ComplexAlignmentGroup):
                        names = []
                        seqs = []
                        regions = []
                        subject = []
                        for align in group.alignments:
                            subject.append(align.subject_key)
                            rec = self.container_factory.seqdb[align.subject_key]
                            seqs.append(
                                str(
                                    rec[
                                        align.subject_region.a : align.subject_region.b
                                    ].seq
                                )
                            )
                            regions.append(
                                (align.subject_region.a, align.subject_region.b)
                            )
                            names.append(rec.name)
                        subject_rec_name = ", ".join(names)
                        subject_seq = ", ".join(seqs)
                        subject_region = regions[:]
                        subject = ",".join(subject)

                    fragments.append(
                        {
                            "query": qk,
                            "query_name": record.name,
                            "query_region": (
                                group.query_region.a,
                                group.query_region.b,
                            ),
                            "subject": subject,
                            "subject_name": subject_rec_name,
                            "subject_region": subject_region,
                            "fragment_length": len(group.query_region),
                            "fragment_seq": subject_seq,
                            "cost": cost,
                            "type": edata["type"],
                        }
                    )

                    # TODO: design overhangs (how long?)
                    # if n1[1]:
                    #     primers.append({
                    #         'query': qk,
                    #         'query_name': record.name,
                    #         'query_region': (align.query_region.a, align.query_region.b),
                    #         'subject': sk,
                    #         'subject_name': subject_rec.name,
                    #         'subject_region': (align.subject_region.a, align.subject_region.a + 20),
                    #         'anneal_seq': str(subject_rec[align.subject_region.a:align.subject_region.a + 20].seq),
                    #         'overhang_seq': '?',
                    #         'cost': '?',
                    #         'type': 'PRIMER'
                    #     })
                    # if n2[1]:
                    #     primers.append({
                    #         'query': qk,
                    #         'query_name': record.name,
                    #         'query_region': (align.query_region.a, align.query_region.b),
                    #         'subject': sk,
                    #         'subject_name': subject_rec.name,
                    #         'subject_region': (align.subject_region.b - 20, align.subject_region.b),
                    #         'fragment_length': 0,
                    #         'anneal_seq': str(subject_rec[align.subject_region.b-20:align.subject_region.b].reverse_complement().seq),
                    #         'overhang_seq': '?',
                    #         'cost': '?',
                    #         'type': 'PRIMER'
                    #     })

                else:
                    B = n1[0]
                    A = n2[0]
                    span = Span(
                        B, A, len(record), cyclic=is_circular(record), allow_wrap=True
                    )

                    # TODO: extending the gene synthesis
                    if not n1[1]:
                        span.b = span.b - 20
                    if not n2[1]:
                        span.a = span.a + 20

                    ranges = span.ranges()
                    frag_seq = record[ranges[0][0] : ranges[0][1]]
                    for r in ranges[1:]:
                        frag_seq += record[r[0] : r[1]]

                    fragments.append(
                        {
                            "query": qk,
                            "query_name": record.name,
                            "query_region": (B, A),
                            "subject": None,
                            "subject_name": "SYNTHESIS",
                            "subject_region": None,
                            "fragment_length": len(span),
                            "fragment_seq": str(frag_seq.seq),
                            "cost": cost,
                            "type": edata["type"],
                        }
                    )
        return pd.DataFrame(fragments), pd.DataFrame(primers)

    # TODO: n_paths to class attribute
    def optimize(self, n_paths=20) -> Dict[str, List[List[AssemblyNode]]]:
        """Finds the optimal paths for each query in the design."""
        results = {}
        for query_key, graph in self.logger.tqdm(
            self.graphs.items(), "INFO", desc="optimizing graphs"
        ):
            container = self.containers[query_key]
            query = container.seqdb[query_key]
            result = DesignResult(container=container, query_key=query_key, graph=graph)
            results[query_key] = result

            cyclic = is_circular(query)
            self.logger.info("Optimizing {}".format(query_key))
            paths = self._collect_optimized_paths(
                graph, len(query), cyclic, n_paths=n_paths
            )
            if not paths:
                query_rec = self.blast_factory.db.records[query_key]
                self.logger.error(
                    "\n\tThere were no solutions found for design '{}' ({}).\n\tThis sequence may"
                    " be better synthesized. Use a tool such as JBEI's BOOST.".format(
                        query_rec.name, query_key
                    )
                )
            result.add_assemblies(paths)
        return results

    def _collect_cycle_endpoints(self, graph: nx.DiGraph, length: int):
        """
        Use the floyd-warshall algorithm to compute the shortest path endpoints.

        :param graph: the networkx graph
        :param length: the size of the query sequence
        :return:
        """
        nodelist, nkeys = sort_with_keys(list(graph.nodes()), key=lambda x: x[0])
        node_to_i = {v: i for i, v in enumerate(nodelist)}
        weight_matrix = np.array(
            nx.floyd_warshall_numpy(graph, nodelist=nodelist, weight="weight")
        )
        endpoints = []

        def bisect_iterator(nodelist, nkeys):
            _i = bisect.bisect_right(nkeys, length)
            for i, A in enumerate(nodelist[:_i]):
                _j = bisect.bisect_left(nkeys, A.index + length)
                for B in nodelist[_j:]:
                    if B.type == "B":
                        j = node_to_i[B]
                        yield i, j, A, B

        pair_iterator = bisect_iterator(nodelist, nkeys)
        for i, j, A, B in pair_iterator:

            # TODO: must include final edge
            a = weight_matrix[i, j]
            b = weight_matrix[j, i]
            if a != np.inf and b != np.inf:
                x = ((A, B), (a, b), a + b)
                endpoints.append(x)

        endpoints = sorted(endpoints, key=lambda x: (x[-1], x[0]))
        return endpoints

    def _nodes_to_fullpaths(
        self,
        graph: nx.DiGraph,
        cycle_endpoints: Tuple[Tuple, Tuple, float],
        cyclic: bool,
        n_paths=None,
    ) -> List[List[Tuple]]:
        """
        Recover full paths from  cycle endpoints.

        :param graph:
        :param cycle_endpoints:
        :param n_paths:
        :return:
        """
        unique_cyclic_paths = []
        for c in cycle_endpoints:
            if n_paths is not None and len(unique_cyclic_paths) >= n_paths:
                break
            path = multipoint_shortest_path(
                graph, c[0], weight_key="weight", cyclic=cyclic
            )
            if path not in unique_cyclic_paths:
                unique_cyclic_paths.append(path)
        return unique_cyclic_paths

    def _collect_optimized_paths(
        self, graph: nx.DiGraph, length: int, cyclic: bool, n_paths=20
    ):
        """
        Collect minimum cycles or linear paths from a graph.

        :param graph: the networkx graph representing the assembly graph
        :param length: length of the query
        :param cyclic: whether to search for a cyclic assembly
        :param n_paths: maximum number of paths to return
        :return:
        """
        if cyclic:
            nodes = self._collect_cycle_endpoints(graph, length=length)
        else:
            raise NotImplementedError("Linear assemblies are not yet implemented.")
        paths = self._nodes_to_fullpaths(graph, nodes, cyclic=cyclic, n_paths=n_paths)
        self._check_paths(paths)
        return paths

    def _check_paths(self, paths):
        """
        Validates a path to check for duplicate nodes.
        :param paths:
        :return:
        """
        invalid_paths = []
        for path in paths:
            lastseen = path[0][2]
            for p in path[1:]:
                if p[2] == lastseen:
                    invalid_paths.append(path)
                    break
                lastseen = p[2]
        if invalid_paths:
            raise DasiDesignException(
                "There are {} invalid paths:\n{}\n...{} more".format(
                    len(invalid_paths),
                    "\n".join([str(x) for x in invalid_paths[:5]]),
                    max(len(invalid_paths) - 5, 0),
                )
            )


class LibraryDesign(Design):
    """
    Design class for producing assemblies for libraries.
    """

    def __init__(self, span_cost=None):
        super().__init__(span_cost)
        self.shared_alignments = []
        self._edges = []

    # @staticmethod
    # def _get_repeats_from_results(results):
    #     repeats = []
    #     for r in results:
    #         qk = r['query']['origin_key']
    #         sk = r['subject']['origin_key']
    #         if qk == sk:
    #             repeats.append((qk, r['query']['start'], r['query']['end']))
    #     return repeats

    def _get_iter_repeats(self, alignments: List[Alignment]):
        """
        Return repeat regions of alignments
        :param alignments:
        :return:
        """
        for align in alignments:
            qk = align.query_key
            sk = align.subject_key
            if qk == sk:
                yield (qk, align.query_region.a, align.query_region.b)

    def _share_query_blast(self):
        """
        Find and use shared fragments across queries.

        :return:
        """
        self.logger.info("=== Expanding shared library fragments ===")

        blast = self.blast_factory(self.QUERIES, self.QUERIES)
        blast.update_config(BLAST_PENALTY_CONFIG)
        blast.quick_blastn()

        results = blast.get_perfect()

        self.logger.info(
            "Found {} shared alignments between the queries".format(len(results))
        )
        self.shared_alignments = results

        self.container_factory.seqdb.update(blast.seq_db.records)
        self.container_factory.load_blast_json(results, Constants.SHARED_FRAGMENT)

        # TODO: expand the normal fragments with the shared fragments
        for query_key, container in self.container_factory.containers().items():
            # expand the share fragments using their own endpoints
            original_shared_fragments = container.get_groups_by_types(
                Constants.SHARED_FRAGMENT
            )
            new_shared_fragments = container.expand_overlaps(
                original_shared_fragments, Constants.SHARED_FRAGMENT
            )

            self.logger.info(
                "{}: Expanded {} shared from original {} shared fragments".format(
                    query_key, len(new_shared_fragments), len(original_shared_fragments)
                )
            )

            # expand the existing fragments with endpoints from the share alignments

            # TODO: what if there is no template for shared fragment?
            # TODO: shared fragment has to be contained wholly in another fragment
            new_alignments = container.expand_overlaps(
                container.get_groups_by_types(
                    [
                        Constants.FRAGMENT,
                        Constants.PCR_PRODUCT,
                        Constants.SHARED_FRAGMENT,
                    ]
                ),
                Constants.PCR_PRODUCT,
            )
            self.logger.info(
                "{}: Expanded {} using {} and found {} new alignments.".format(
                    query_key,
                    Constants.PCR_PRODUCT,
                    Constants.SHARED_FRAGMENT,
                    len(new_alignments),
                )
            )
            # grab the pcr products and expand primer pairs (again)
            templates = container.get_groups_by_types(Constants.PCR_PRODUCT)
            new_primer_pairs = container.expand_primer_pairs(templates)
            self.logger.info(
                "{}: Expanded {} {} using {}".format(
                    query_key,
                    len(new_primer_pairs),
                    "PRODUCTS_WITH_PRIMERS",
                    Constants.SHARED_FRAGMENT,
                )
            )

        repeats = []
        for query_key, container in self.container_factory.containers().items():
            # get all shared fragments
            alignments = container.get_alignments_by_types(Constants.SHARED_FRAGMENT)
            self.logger.info(
                "{} shared fragments for {}".format(len(alignments), query_key)
            )
            # add to list of possible repeats
            repeats += list(self._get_iter_repeats(alignments))
        self.repeats = repeats

    def compile_library(self):
        """Compile the materials list into assembly graphs."""
        self.graphs = {}
        self._blast()
        self._share_query_blast()
        self.assemble_graphs()

    def optimize_library(self):
        """Optimize the assembly graph for library assembly."""
        raise NotImplementedError
