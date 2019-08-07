"""Primer and synthesis design"""

from shoestring import AlignmentContainer, Constants, AssemblyGraphBuilder
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
        self.factory = BioBlastFactory()
        self.logger = logger(self)
        self.G = None
        self.span_cost = span_cost
        self.container = None

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
        self.factory.add_records(primers, self.PRIMERS)

    def add_templates(self, templates: List[SeqRecord]):
        self.logger.info("Adding templates")
        self.factory.add_records(templates, self.TEMPLATES)

    def add_queries(self, queries: List[SeqRecord]):
        self.logger.info("Adding queries")
        self.factory.add_records(queries, self.QUERIES)

    def _blast(self):
        self.logger.info("Compiling assembly graph")

        blast = self.factory("templates", "queries")
        primer_blast = self.factory("primers", "queries")

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

        self.container = AlignmentContainer(seqdb)

        # load the results (potential PCR Products)
        self.container.load_blast_json(results, Constants.PCR_PRODUCT)

        # load the primer results
        self.container.load_blast_json(primer_results, Constants.PRIMER)

    def compile(self):
        self._blast()

        self.container.expand(expand_overlaps=True, expand_primers=True)

        # group by query_regions
        groups = self.container.alignment_groups

        self.logger.info("Number of types: {}".format(len(self.container.groups_by_type)))
        self.logger.info("Number of groups: {}".format(len(groups)))

        # build assembly graph
        graph_builder = AssemblyGraphBuilder(self.container, span_cost=self.span_cost)
        G = graph_builder.build_assembly_graph()

        self.logger.info("=== Assembly Graph ===")
        self.logger.info(nx.info(G))
        assert G.number_of_edges()
        self.G = G

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
        self.logger.info("Optimizing...")

        # shortest path matrix
        nodelist = list(self.G.nodes())
        weight_matrix = np.array(nx.floyd_warshall_numpy(self.G, nodelist=nodelist, weight='weight'))

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
            path1 = nx.shortest_path(self.G, c[0], c[1], weight='weight')
            path2 = nx.shortest_path(self.G, c[1], c[0], weight='weight')
            path = path1 + path2[1:]
            for n1, n2 in pairwise(path):
                edata = self.G[n1][n2]
                print('{} > {} Weight={} name={} span={}'.format(n1, n2, edata['weight'], edata['name'], edata['span_length']))
            print()