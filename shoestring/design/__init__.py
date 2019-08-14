"""Primer and synthesis design"""

from shoestring import AlignmentContainerFactory, Constants, AssemblyGraphBuilder
from shoestring.utils import perfect_subject
import networkx as nx
from shoestring import BioBlastFactory
from shoestring.log import logger
from typing import List
from Bio.SeqRecord import SeqRecord
import numpy as np
from more_itertools import pairwise
from abc import ABC, abstractmethod
from pyblast.utils import Span
import pandas as pd


class DNADesign(ABC):

    @abstractmethod
    def design(self, span, seq1, seq2, query_seq):
        pass


class PrimerDesign(object):

    # TODO: how much to expand?
    # TODO: what kind of constraints?
    # TODO: re-evaluate costs after refinement
    def design(self, span, seq1, seq2, query_seq, expand_left=True, expand_right=True):
        pass


class SynthesisDesign(object):

    def design(self, span, seq1, seq2, query_seq):
        pass


class Design(object):

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"

    def __init__(self, span_cost=None):
        self.blast_factory = BioBlastFactory()
        self.logger = logger(self)

        # graph by query_key
        self.graphs = {}
        self.span_cost = span_cost
        self.container_factory = None

    def add_materials(
        self,
        primers: List[SeqRecord],
        templates: List[SeqRecord],
        queries: List[SeqRecord],
    ):
        self.add_primers(primers)
        self.add_templates(templates)
        self.add_queries(queries)

    def add_primers(self, primers: List[SeqRecord]):
        self.logger.info("Adding primers")
        self.blast_factory.add_records(primers, self.PRIMERS)

    def add_templates(self, templates: List[SeqRecord]):
        self.logger.info("Adding templates")
        self.blast_factory.add_records(templates, self.TEMPLATES)

    def add_queries(self, queries: List[SeqRecord]):
        self.logger.info("Adding queries")
        self.blast_factory.add_records(queries, self.QUERIES)

    def _blast(self):
        self.logger.info("Compiling assembly graph")

        blast = self.blast_factory("templates", "queries")
        primer_blast = self.blast_factory("primers", "queries")

        blast.quick_blastn()
        primer_blast.quick_blastn_short()

        results = blast.get_perfect()
        primer_results = primer_blast.get_perfect()

        primer_results = [p for p in primer_results if perfect_subject(p["subject"])]
        self.logger.info("Number of perfect primers: {}".format(len(primer_results)))
        # primer_results = [p for p in primer_results if p['subject']['start'] == 1]

        # combine the sequence databases (for now)
        seqdb = {}
        seqdb.update(blast.seq_db.records)
        seqdb.update(primer_blast.seq_db.records)

        self.container_factory = AlignmentContainerFactory(seqdb)
        self.container_factory.load_blast_json(results, Constants.PCR_PRODUCT)
        self.container_factory.load_blast_json(primer_results, Constants.PRIMER)

    def container_list(self):
        return list(self.container_factory.containers().values())

    def compile(self):
        """Compile materials to assembly graph"""
        self._blast()

        for query_key, container in self.logger.tqdm(self.container_factory.containers().items(), "INFO", desc='compiling all containers'):
            container.expand(expand_overlaps=True, expand_primers=True)

            # group by query_regions
            groups = container.alignment_groups

            self.logger.info("Number of types: {}".format(len(container.groups_by_type)))
            self.logger.info("Number of groups: {}".format(len(groups)))

            # build assembly graph
            graph_builder = AssemblyGraphBuilder(container, span_cost=self.span_cost)
            G = graph_builder.build_assembly_graph()

            self.logger.info("=== Assembly Graph ===")
            self.logger.info(nx.info(G))
            assert G.number_of_edges()
            self.graphs[query_key] = G

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

    def path_to_df(self, paths_dict):
        def find(a, b, alignments):
            for align in alignments:
                if a == align.query_region.a and b == align.query_region.b:
                    yield align

        rows = []

        for qk, paths in paths_dict.items():
            G = self.graphs[qk]
            alignments = self.container_factory.alignments[qk]
            record = self.container_factory.seqdb[qk]
            path = paths[0]

            for n1, n2 in pairwise(path):
                edata = G[n1][n2]
                cost = edata['weight']
                print(edata)
                if n1[-1] == 'A' and n2[-1] == 'B':
                    A = n1[0]
                    B = n2[0]
                    align = list(find(A, B, alignments))[0]
                    sk = align.subject_key
                    subject_rec = self.container_factory.seqdb[sk]
                    subject_seq = str(subject_rec[align.subject_region.a:align.subject_region.b].seq)

                    rows.append({
                        'query': qk,
                        'query_name': record.name,
                        'query_region': (align.query_region.a, align.query_region.b),
                        'subject': sk,
                        'subject_name': subject_rec.name,
                        'subject_region': (align.subject_region.a, align.subject_region.b),
                        'fragment_length': len(align.subject_region),
                        'fragment_seq': subject_seq,
                        'cost': cost,
                        'type': edata['type']
                    })
                else:
                    B = n1[0]
                    A = n2[0]
                    span = Span(B, A, len(record), cyclic=is_circular(record), allow_wrap=True)
                    ranges = span.ranges()
                    frag_seq = record[ranges[0][0]:ranges[0][1]]
                    for r in ranges[1:]:
                        frag_seq += record[r[0]:r[1]]

                    rows.append({
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
        pd.DataFrame(rows)

    def design(self):
        path_dict = self.optimize()
        df = self.path_to_df(path_dict)
        return df

    def optimize(self, verbose=False):
        query_key_to_path = {}
        for query_key, G in self.logger.tqdm(self.graphs.items(), "INFO", desc='optimizing graphs'):
            self.logger.info("Optimizing {}".format(query_key))
            paths = self._optimize_graph(G)
            if verbose:
                for path in paths:
                    for n1, n2 in pairwise(path):
                        edata = G[n1][n2]
                        print('{} > {} Weight={} name={} span={} type={}'.format(n1, n2, edata['weight'], edata['name'],                                                      edata['span_length'], edata['type']))
                    print()
            query_key_to_path[query_key] = paths
        return query_key_to_path

    def path_to_design(self, graph, query_key):
        fragments = []
        for n1, n2, edata in graph.edges(data=True):
            if n1[0] == 'A' and n2[0] == 'B':
                pass

    def _optimize_graph(self, graph):

        # shortest path matrix
        nodelist = list(graph.nodes())
        weight_matrix = np.array(nx.floyd_warshall_numpy(graph, nodelist=nodelist, weight='weight'))

        # shortest cycles (estimated)
        cycles = []
        paths = []
        for i, _ in enumerate(weight_matrix):
            for j, _ in enumerate(weight_matrix[i]):
                a = weight_matrix[i, j]
                b = weight_matrix[j, i]
                if i == j:
                    continue

                anode = nodelist[i]
                bnode = nodelist[j]
                if a != np.inf:
                    paths.append((anode, bnode, a))
                if b != np.inf:
                    paths.append((bnode, anode, b))
                if a != np.inf and b != np.inf:
                    cycles.append((anode, bnode, a, b, a + b))

        cycles = sorted(cycles, key=lambda c: c[-1])

        self.logger.info("Cycles: {}".format(len(cycles)))
        self.logger.info("Paths: {}".format(len(paths)))

        # print cycles
        paths = []
        for c in cycles[:20]:
            print(c)
            path1 = nx.shortest_path(graph, c[0], c[1], weight='weight')
            path2 = nx.shortest_path(graph, c[1], c[0], weight='weight')
            path = path1 + path2[1:]
            paths.append(path)
        return paths