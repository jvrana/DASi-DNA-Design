"""Primer and synthesis design

.. module:: design

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    assembly
"""
from __future__ import annotations

import bisect
from collections.abc import Iterable
from copy import deepcopy
from itertools import zip_longest
from typing import List, Tuple, Dict

import networkx as nx
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from pyblast import BioBlastFactory
from pyblast.utils import Span, is_circular
from multiprocessing import Pool
from dasi.alignments import (
    Alignment,
    AlignmentContainerFactory,
    AlignmentContainer,
    AlignmentGroup,
    ComplexAlignmentGroup,
)
from dasi.constants import Constants
from dasi.design.graph_builder import AssemblyNode
from dasi.log import logger
from dasi.utils import perfect_subject, sort_cycle
from .design_algorithms import assemble_graph, optimize_graph, multiprocessing_assemble_graph, \
    multiprocessing_optimize_graph
from .graph_builder import AssemblyGraphBuilder

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
        self.groups = container.groups()
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
                    "cost": np.inf,
                    "weight": np.inf,
                    "type": "missing",
                    "span": np.inf,
                    "name": "missing",
                }

            # TODO: fix query_region (overlaps are backwards)
            query_region = self.container.alignments[0].query_region.new(
                n1.index, n2.index, allow_wrap=True
            )
            groups = self.container.find_groups_by_pos(
                query_region.a, query_region.b, groups=self.groups
            )
            edata["groups"] = groups
            edata["query_region"] = query_region
            SG.add_edge(
                _resolve(n1, query_region)[0], _resolve(n2, query_region)[0], **edata
            )
        return SG

    @property
    def cyclic(self):
        return is_circular(self.query)

    # TODO: this cost is no longer true...
    def cost(self):
        total = 0
        for _, _, edata in self.edges():
            if "cost" not in edata:
                x = 1
            total += edata["cost"]
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
                    "cost": edata["cost"],
                    "span": edata["span"],
                    "type": edata["type"],
                    "name": edata["name"],
                    "efficiency": edata.get("efficiency", np.nan),
                }
            )

        df = pd.DataFrame(rows)
        return df

    def __eq__(self, other: Assembly) -> bool:
        return self.edit_distance(other) == 0

    def __iter__(self):
        for n in self.nodes(data=False):
            yield n

def run_blast(args):
    if args is None or args[0] is None:
        return []
    blast, func_str = args
    func = getattr(blast, func_str)
    func()
    return blast.get_perfect()


class FakePool(object):

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return FakePool()

    def __exit__(self, a, b, c):
        pass

    def map(self, func, args):
        return [func(arg) for arg in args]


class Design(object):
    """
    Design class that returns optimal assemblies from a set of materials.
    """

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"
    FRAGMENTS = "fragments"
    DEFAULT_N_THREADS = 1

    def __init__(self, span_cost=None, seqdb=None, n_threads=None):
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
        self.n_threads = n_threads or self.DEFAULT_N_THREADS

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
        """Add primer sequences to materials list"""
        self.logger.info("Adding primers")
        self.blast_factory.add_records(primers, self.PRIMERS)

    def add_templates(self, templates: List[SeqRecord]):
        """Add template sequences to materials list"""
        self.logger.info("Adding templates")
        self.blast_factory.add_records(templates, self.TEMPLATES)

    def add_queries(self, queries: List[SeqRecord]):
        """Add goal/query sequences to materials list"""
        self.logger.info("Adding queries")
        self.blast_factory.add_records(queries, self.QUERIES)

    def add_fragments(self, fragments: List[SeqRecord]):
        """Add fragment sequences to materials list"""
        self.logger.info("Adding fragments")
        self.blast_factory.add_records(fragments, self.FRAGMENTS)
        # self.blast_factory.add_records(fragments, self.TEMPLATES)

    @classmethod
    def filter_linear_records(cls, records):
        """Return only linear records"""
        return [r for r in records if not is_circular(r)]

    @classmethod
    def filter_perfect_subject(cls, results):
        """return only results whose subject is 100% aligned to query"""
        return [r for r in results if perfect_subject(r["subject"])]

    # # TODO: do a single blast and sort results based on record keys
    # def _blast(self):
    #     """Preform blast of materials against queries."""
    #     self.logger.info("Compiling assembly graph")
    #
    #     with Pool(processes=3) as pool:
    #         template_blast = self.blast_factory(self.TEMPLATES, self.QUERIES)
    #         template_blast.update_config(BLAST_PENALTY_CONFIG)
    #         self.container_factory.seqdb.update(template_blast.seq_db.records)
    #
    #         if self.blast_factory.record_groups[self.FRAGMENTS]:
    #             fragment_blast = self.blast_factory(self.FRAGMENTS, self.QUERIES)
    #             fragment_blast.update_config(BLAST_PENALTY_CONFIG)
    #             self.container_factory.seqdb.update(fragment_blast.seq_db.records)
    #         else:
    #             fragment_blast = None
    #
    #         if self.blast_factory.record_groups[self.PRIMERS]:
    #             primer_blast = self.blast_factory(self.PRIMERS, self.QUERIES)
    #             primer_blast.update_config(BLAST_PENALTY_CONFIG)
    #             self.container_factory.seqdb.update(primer_blast.seq_db.records)
    #         else:
    #             primer_blast = None
    #
    #         template_results, fragment_results, primer_results = pool.map(run_blast, (
    #             (template_blast, 'quick_blastn'),
    #             (fragment_blast, 'quick_blastn'),
    #             (primer_blast, 'quick_blastn_short')
    #         ))
    #
    #         fragment_results = self.filter_perfect_subject(fragment_results)
    #         primer_results = self.filter_perfect_subject(primer_results)
    #
    #         self.container_factory.load_blast_json(fragment_results, Constants.FRAGMENT)
    #         self.container_factory.load_blast_json(template_results, Constants.PCR_PRODUCT)
    #         self.container_factory.load_blast_json(primer_results, Constants.PRIMER)

    # non-threaded blast
    def _blast(self):
        # align templates
        template_blast = self.blast_factory(self.TEMPLATES, self.QUERIES)
        template_blast.update_config(BLAST_PENALTY_CONFIG)
        template_blast.quick_blastn()
        template_results = template_blast.get_perfect()
        self.template_results = template_results
        self.logger.info("Number of template matches: {}".format(len(template_results)))
        self.container_factory.seqdb.update(template_blast.seq_db.records)

        # align fragments
        if self.blast_factory.record_groups[self.FRAGMENTS]:
            fragment_blast = self.blast_factory(self.FRAGMENTS, self.QUERIES)
            fragment_blast.update_config(BLAST_PENALTY_CONFIG)
            fragment_blast.quick_blastn()
            fragment_results = self.filter_perfect_subject(fragment_blast.get_perfect())
            self.container_factory.seqdb.update(fragment_blast.seq_db.records)
            self.logger.info(
                "Number of perfect fragment matches: {}".format(len(fragment_results))
            )
        else:
            fragment_results = []
        self.fragment_results = fragment_results

        # align primers
        if self.blast_factory.record_groups[self.PRIMERS]:
            primer_blast = self.blast_factory(self.PRIMERS, self.QUERIES)
            primer_blast.update_config(BLAST_PENALTY_CONFIG)
            primer_blast.quick_blastn_short()
            primer_results = self.filter_perfect_subject(primer_blast.get_perfect())
            self.container_factory.seqdb.update(primer_blast.seq_db.records)
            self.logger.info(
                "Number of perfect primers: {}".format(len(primer_results))
            )
        else:
            primer_results = []
        self.primer_results = primer_results

        self.container_factory.load_blast_json(fragment_results, Constants.FRAGMENT)
        self.container_factory.load_blast_json(template_results, Constants.PCR_PRODUCT)
        self.container_factory.load_blast_json(primer_results, Constants.PRIMER)

    @property
    def containers(self):
        """Iterable of alignment containers in this design."""
        return self.container_factory.containers()

    def container_list(self):
        """List of alignment containers in this design."""
        return list(self.container_factory.containers().values())

    def query_keys(self):
        """List of query keys in this design."""
        return list(self.container_factory.containers())

    def assemble_graphs(self, n_cores=None):
        n_cores = n_cores or self.n_threads
        if n_cores > 1:
            with self.logger.timeit("INFO",
                                    "assembling graphs (n_graphs={}, threads={})".format(len(self.container_list()),
                                                                                         n_cores)):
                return self._assemble_graphs_with_threads(n_cores)
        return self._assemble_graphs_without_threads()

    def _assemble_graphs_without_threads(self):
        """Assemble all assembly graphs for all queries in this design."""
        for query_key, container in self.logger.tqdm(
                self.container_factory.containers().items(),
                "INFO",
                desc="assembling graphs (threads=1)",
        ):
            self.graphs[query_key] = assemble_graph(container, self.span_cost)

    def _assemble_graphs_with_threads(self, n_cores=None):
        query_keys, containers = zip(*self.container_factory.containers().items())
        graphs = multiprocessing_assemble_graph(containers, self.span_cost, n_cores=n_cores)
        for qk, g in zip(query_keys, graphs):
            self.graphs[qk] = g

    def compile(self, n_cores=None):
        """Compile materials to assembly graph"""
        self.results = {}
        with self.logger.timeit("INFO", "running blast"):
            self._blast()
        self.assemble_graphs(n_cores=n_cores)

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
                            align.subject_region.a: align.subject_region.b
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
                                    align.subject_region.a: align.subject_region.b
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
                    frag_seq = record[ranges[0][0]: ranges[0][1]]
                    for r in ranges[1:]:
                        frag_seq += record[r[0]: r[1]]

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

    def optimize(self, n_paths=3, n_cores=None):
        n_cores = n_cores or self.n_threads
        if n_cores > 1:
            with self.logger.timeit("INFO",
                                    "optimizing graphs (n_graphs={}, threads={})".format(len(self.graphs), n_cores)):
                return self._optimize_with_threads(n_paths, n_cores)
        return self._optimize_without_threads(n_paths)

    # TODO: n_paths to class attribute
    def _optimize_without_threads(self, n_paths=5) -> Dict[str, List[List[AssemblyNode]]]:
        """Finds the optimal paths for each query in the design."""
        results_dict = {}
        for query_key, graph, query_length, cyclic, result in self.logger.tqdm(
                self._collect_optimize_args(self.graphs),
                "INFO",
                desc="optimizing graphs (n_graphs={}, threads=1)".format(len(self.graphs))
        ):
            container = self.containers[query_key]
            query = container.seqdb[query_key]
            cyclic = is_circular(query)
            results_dict[query_key] = result
            paths = optimize_graph(graph, len(query), cyclic, n_paths)
            if not paths:
                query_rec = self.blast_factory.db.records[query_key]
                self.logger.error(
                    "\n\tThere were no solutions found for design '{}' ({}).\n\tThis sequence may"
                    " be better synthesized. Use a tool such as JBEI's BOOST.".format(
                        query_rec.name, query_key
                    )
                )
            result.add_assemblies(paths)
        return results_dict

    def _optimize_with_threads(self, n_paths=5, n_cores=10):
        results_dict = {}
        query_keys, graphs, query_lengths, cyclics, results = zip(*list(self._collect_optimize_args(self.graphs)))

        list_of_paths = multiprocessing_optimize_graph(graphs=graphs, query_lengths=query_lengths, cyclics=cyclics,
                                                       n_paths=n_paths, n_cores=n_cores)
        for qk, paths, result in zip(query_keys, list_of_paths, results):
            result.add_assemblies(paths)
            results_dict[qk] = result
        return results_dict

    def _collect_optimize_args(self, graphs):
        for query_key, graph in self.logger.tqdm(
                graphs.items(), "INFO", desc="optimizing graphs"
        ):
            container = self.containers[query_key]
            query = container.seqdb[query_key]
            cyclic = is_circular(query)
            result = DesignResult(container=container, query_key=query_key, graph=graph)
            yield query_key, graph, len(query), cyclic, result

    # TODO: make this a method outside of class scope for multithreading.
    #       In order to do this, all of the methods will need to be scoped
    #       outside of the class, probably in its own file as 'design_algorithms.py'


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
