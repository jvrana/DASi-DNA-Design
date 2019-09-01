"""Primer and synthesis design

.. module:: design

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    assembly
"""

from dasi.alignments import Alignment, AlignmentContainerFactory, AlignmentGroup, ComplexAlignmentGroup
from dasi.constants import Constants
from .assembly import AssemblyGraphBuilder
from dasi.utils import perfect_subject, multipoint_shortest_path, sort_with_keys
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

BLAST_PENALTY_CONFIG = {
    'gapopen': 3,
    'gapextend': 3,
    'reward': 1,
    'penalty': -5
}


class DesignBase(object):
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
        self.graphs = {}
        self.span_cost = span_cost
        self.container_factory = AlignmentContainerFactory(self.seqdb)

    @property
    def seqdb(self):
        return self._seqdb


class DesignResult(object):
    """
    Should take in a path, graph, container, seqdb to produce relevant information
    """

    def __init__(self, nodes):
        self.nodes = nodes


class Design(DesignBase):
    """
    Design class that returns optimal assemblies from a set of materials.
    """

    def add_materials(
            self,
            primers: List[SeqRecord],
            templates: List[SeqRecord],
            queries: List[SeqRecord],
            fragments=None
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
        self.logger.info("Number of perfect fragment matches: {}".format(len(fragment_results)))

        # align primers
        if self.blast_factory.record_groups[self.PRIMERS]:
            primer_blast = self.blast_factory(self.PRIMERS, self.QUERIES)
            primer_blast.update_config(BLAST_PENALTY_CONFIG)
            primer_blast.quick_blastn_short()
            primer_results = primer_blast.get_perfect()
            primer_results = self.filter_perfect_subject(primer_results)
            self.container_factory.seqdb.update(primer_blast.seq_db.records)
            self.logger.info("Number of perfect primers: {}".format(len(primer_results)))
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

    def assemble_graphs(self):
        for query_key, container in self.logger.tqdm(self.container_factory.containers().items(), "INFO",
                                                     desc='compiling all containers'):
            container.expand(expand_overlaps=True, expand_primers=True)

            # group by query_regions
            groups = container.groups()

            self.logger.info("Number of types: {}".format(len(container.groups_by_type)))
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
        self.graphs = {}
        self._blast()
        self.assemble_graphs()

    # def plot_matrix(self, matrix):
    ## plot matrix
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
            sub = record[ranges[0][0]:ranges[0][1]]
            for r in ranges[1:]:
                sub += record[r[0]:r[1]]
            sub.annotations = record.annotations
            return sub

        alignments = self.container_factory.alignments[query_key]
        align = list(self._find_iter_alignment(a, b, alignments))[0]
        subject_key = align.subject_key
        subject_rec = self.container_factory.seqdb[subject_key]
        query_rec = self.container_factory.seqdb[query_key]

        subject_seq = sub_record(subject_rec, align.subject_region)

        fragment_info = {
            'query_id': query_key,
            'query_name': query_rec.name,
            'query_region': (align.query_region.a, align.query_region.b),
            'subject_id': subject_key,
            'subject_name': subject_rec.name,
            'subject_region': (align.subject_region.a, align.subject_region.b),
            'fragment_length': len(align.subject_region),
            'fragment_seq': subject_seq,
            'cost': cost,
            'type': fragment_type,
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
                cost = edata['weight']
                if n1[2] == 'A' and n2[2] == 'B':
                    A = n1[0]
                    B = n2[0]
                    group = container.find_groups_by_pos(A, B)[0]

                    if isinstance(group, AlignmentGroup):
                        align = group.alignments[0]
                        subject = align.subject_key
                        subject_rec = self.container_factory.seqdb[align.subject_key]
                        subject_rec_name = subject_rec.name
                        subject_seq = str(subject_rec[align.subject_region.a:align.subject_region.b].seq)
                        subject_region = (align.subject_region.a, align.subject_region.b)
                    elif isinstance(group, ComplexAlignmentGroup):
                        names = []
                        seqs = []
                        regions = []
                        subject = []
                        for align in group.alignments:
                            subject.append(align.subject_key)
                            rec = self.container_factory.seqdb[align.subject_key]
                            seqs.append(str(rec[align.subject_region.a:align.subject_region.b].seq))
                            regions.append((align.subject_region.a, align.subject_region.b))
                            names.append(rec.name)
                        subject_rec_name = ', '.join(names)
                        subject_seq = ', '.join(seqs)
                        subject_region = regions[:]
                        subject = ','.join(subject)

                    fragments.append({
                        'query': qk,
                        'query_name': record.name,
                        'query_region': (group.query_region.a, group.query_region.b),
                        'subject': subject,
                        'subject_name': subject_rec_name,
                        'subject_region': subject_region,
                        'fragment_length': len(group.query_region),
                        'fragment_seq': subject_seq,
                        'cost': cost,
                        'type': edata['type']
                    })

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
                    span = Span(B, A, len(record), cyclic=is_circular(record), allow_wrap=True)

                    # TODO: extending the gene synthesis
                    if not n1[1]:
                        span.b = span.b - 20
                    if not n2[1]:
                        span.a = span.a + 20

                    ranges = span.ranges()
                    frag_seq = record[ranges[0][0]:ranges[0][1]]
                    for r in ranges[1:]:
                        frag_seq += record[r[0]:r[1]]

                    fragments.append({
                        'query': qk,
                        'query_name': record.name,
                        'query_region': (B, A),
                        'subject': None,
                        'subject_name': 'SYNTHESIS',
                        'subject_region': None,
                        'fragment_length': len(span),
                        'fragment_seq': str(frag_seq.seq),
                        'cost': cost,
                        'type': edata['type']
                    })
        return pd.DataFrame(fragments), pd.DataFrame(primers)

    def design(self):
        path_dict = self.optimize()
        df = self.path_to_df(path_dict)
        return df

    # TODO: n_paths to class attribute
    def optimize(self, n_paths=20) -> Dict[str, List[List[Tuple]]]:
        """Finds the optimal paths for each query in the design."""
        query_key_to_path = {}
        for query_key, G in self.logger.tqdm(self.graphs.items(), "INFO", desc='optimizing graphs'):
            container = self.containers[query_key]
            query = container.seqdb[query_key]
            cyclic = is_circular(query)
            self.logger.info("Optimizing {}".format(query_key))
            paths = self._collect_optimized_paths(G, len(query), cyclic, n_paths=n_paths)
            if not paths:
                query_rec = self.blast_factory.db.records[query_key]
                self.logger.error("\n\tThere were no solutions found for design '{}' ({}).\n\tThis sequence may"
                                  " be better synthesized. Use a tool such as JBEI's BOOST.".format(query_rec.name,
                                                                                                    query_key))
            query_key_to_path[query_key] = paths
        return query_key_to_path

    def _collect_cycle_endpoints(self, graph: nx.DiGraph, length: int):
        """
        Use the floyd warshall algorithm to compute the shortest path endpoints.

        :param graph: the networkx graph
        :param length: the size of the query sequence
        :return:
        """
        nodelist, nkeys = sort_with_keys(list(graph.nodes()), key=lambda x: x[0])
        node_to_i = {v: i for i, v in enumerate(nodelist)}
        weight_matrix = np.array(nx.floyd_warshall_numpy(graph, nodelist=nodelist, weight='weight'))
        cycle_endpoints = []

        def bisect_iterator(nodelist, nkeys):
            for i, A in enumerate(nodelist):
                _j = bisect.bisect_left(nkeys, A[0] + length)
                for B in nodelist[_j:]:
                    j = node_to_i[B]
                    yield i, j, A, B

        pair_iterator = bisect_iterator(nodelist, nkeys)
        for i, j, A, B in pair_iterator:
            if B[2] == 'B' and A[2] == 'A':
                a = weight_matrix[i, j]
                b = weight_matrix[j, i]
                if a != np.inf and b != np.inf:
                    x = ((A, B), (a, b), a + b)
                    cycle_endpoints.append(x)

        cycle_endpoints = sorted(cycle_endpoints, key=lambda x: (x[-1], x[0]))
        return cycle_endpoints

    def _nodes_to_fullpaths(self, graph: nx.DiGraph, cycle_endpoints: Tuple[Tuple, Tuple, float], cyclic: bool, n_paths=None) -> \
            List[List[Tuple]]:
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
            path = multipoint_shortest_path(graph, c[0], weight_key='weight', cyclic=cyclic)
            if path not in unique_cyclic_paths:
                unique_cyclic_paths.append(path)
        return unique_cyclic_paths

    def _collect_optimized_paths(self, graph: nx.DiGraph, length: int, cyclic: bool, n_paths=20):
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
            raise DasiDesignException("There are {} invalid paths:\n{}\n...{} more".format(
                len(invalid_paths),
                "\n".join([str(x) for x in invalid_paths[:5]]),
                max(len(invalid_paths) - 5, 0)
            ))


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

        self.logger.info("Found {} shared alignments between the queries".format(len(results)))
        self.shared_alignments = results

        self.container_factory.seqdb.update(blast.seq_db.records)
        self.container_factory.load_blast_json(results, Constants.SHARED_FRAGMENT)

        # TODO: expand the normal fragments with the shared fragments
        for query_key, container in self.container_factory.containers().items():
            # expand the share fragments using their own endpoints
            original_shared_fragments = container.get_groups_by_types(Constants.SHARED_FRAGMENT)
            new_shared_fragments = container.expand_overlaps(original_shared_fragments,
                                                             Constants.SHARED_FRAGMENT)

            self.logger.info("{}: Expanded {} shared from original {} shared fragments".format(
                query_key,
                len(new_shared_fragments),
                len(original_shared_fragments)
            ))

            # expand the existing fragments with endpoints from the share alignments

            # TODO: what if there is no template for shared fragment?
            # TODO: shared fragment has to be contained wholly in another fragment
            new_alignments = container.expand_overlaps(container.get_groups_by_types(
                [Constants.FRAGMENT,
                 Constants.PCR_PRODUCT,
                 Constants.SHARED_FRAGMENT]
            ), Constants.PCR_PRODUCT)
            self.logger.info("{}: Expanded {} using {} and found {} new alignments.".format(
                query_key,
                Constants.PCR_PRODUCT,
                Constants.SHARED_FRAGMENT,
                len(new_alignments)
            ))
            # grab the pcr products and expand primer pairs (again)
            templates = container.get_groups_by_types(
                Constants.PCR_PRODUCT
            )
            new_primer_pairs = container.expand_primer_pairs(templates)
            self.logger.info("{}: Expanded {} {} using {}".format(
                query_key,
                len(new_primer_pairs),
                "PRODUCTS_WITH_PRIMERS",
                Constants.SHARED_FRAGMENT
            ))

        repeats = []
        for query_key, container in self.container_factory.containers().items():
            # get all shared fragments
            alignments = container.get_alignments_by_types(Constants.SHARED_FRAGMENT)
            self.logger.info("{} shared fragments for {}".format(len(alignments), query_key))
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
