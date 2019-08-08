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

    def optimize(self):
        for query_key, G in self.logger.tqdm(self.graphs.items(), "INFO", desc='optimizing graphs'):
            self.logger.info("Optimizing {}".format(query_key))
            self._optimize_graph(G)

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
        for c in cycles[:20]:
            print(c)
            path1 = nx.shortest_path(graph, c[0], c[1], weight='weight')
            path2 = nx.shortest_path(graph, c[1], c[0], weight='weight')
            path = path1 + path2[1:]
            for n1, n2 in pairwise(path):
                edata = graph[n1][n2]
                print('{} > {} Weight={} name={} span={}'.format(n1, n2, edata['weight'], edata['name'], edata['span_length']))
            print()