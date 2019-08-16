"""Primer and synthesis design"""

from shoestring import Alignment, AlignmentContainer, AlignmentContainerFactory, Constants, AssemblyGraphBuilder
from shoestring.utils import perfect_subject
import networkx as nx
from pyblast import BioBlastFactory
from shoestring.log import logger
from typing import List
from Bio.SeqRecord import SeqRecord
import numpy as np
from more_itertools import pairwise
from abc import ABC, abstractmethod
from pyblast.utils import Span, is_circular
import pandas as pd
from typing import Tuple

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


class DesignBase(object):

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"

    def __init__(self, span_cost=None):
        self.blast_factory = BioBlastFactory()
        self.logger = logger(self)

        # graph by query_key
        self.graphs = {}
        self.span_cost = span_cost
        self.container_factory = AlignmentContainerFactory({})


class Design(DesignBase):

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

        self.container_factory.seqdb.update(blast.seq_db.records)
        self.container_factory.seqdb.update(primer_blast.seq_db.records)
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


    def path_to_df(self, paths_dict):
        def find(a, b, alignments):
            for align in alignments:
                if a == align.query_region.a and b == align.query_region.b:
                    yield align

        fragments = []
        primers = []

        for qk, paths in paths_dict.items():
            G = self.graphs[qk]
            alignments = self.container_factory.alignments[qk]
            record = self.container_factory.seqdb[qk]
            path = paths[0]

            for n1, n2 in pairwise(path):
                edata = G[n1][n2]
                cost = edata['weight']
                if n1[-1] == 'A' and n2[-1] == 'B':
                    A = n1[0]
                    B = n2[0]
                    align = list(find(A, B, alignments))[0]
                    sk = align.subject_key
                    subject_rec = self.container_factory.seqdb[sk]
                    subject_seq = str(subject_rec[align.subject_region.a:align.subject_region.b].seq)

                    fragments.append({
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

                    # TODO: design overhangs (how long?)
                    if n1[1]:
                        primers.append({
                            'query': qk,
                            'query_name': record.name,
                            'query_region': (align.query_region.a, align.query_region.b),
                            'subject': sk,
                            'subject_name': subject_rec.name,
                            'subject_region': (align.subject_region.a, align.subject_region.a + 20),
                            'anneal_seq': str(subject_rec[align.subject_region.a:align.subject_region.a + 20].seq),
                            'overhang_seq': '?',
                            'cost': '?',
                            'type': 'PRIMER'
                        })
                    if n2[1]:
                        primers.append({
                            'query': qk,
                            'query_name': record.name,
                            'query_region': (align.query_region.a, align.query_region.b),
                            'subject': sk,
                            'subject_name': subject_rec.name,
                            'subject_region': (align.subject_region.b - 20, align.subject_region.b),
                            'fragment_length': 0,
                            'anneal_seq': str(subject_rec[align.subject_region.b-20:align.subject_region.b].reverse_complement().seq),
                            'overhang_seq': '?',
                            'cost': '?',
                            'type': 'PRIMER'
                        })

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

    def optimize(self, verbose=False):
        query_key_to_path = {}
        for query_key, G in self.logger.tqdm(self.graphs.items(), "INFO", desc='optimizing graphs'):
            self.logger.info("Optimizing {}".format(query_key))
            paths = self._optimize_graph(G)
            if verbose:
                for path in paths:
                    for n1, n2 in pairwise(path):
                        edata = G[n1][n2]
                        print('{} > {} Weight={} name={} span={} type={}'.format(n1, n2, edata['weight'], edata['name']))
            query_key_to_path[query_key] = paths
        return query_key_to_path

    def _optimize_graph(self, graph):

        # shortest path matrix
        nodelist = list(graph.nodes())
        weight_matrix = np.array(nx.floyd_warshall_numpy(graph, nodelist=nodelist, weight='weight'))

        # shortest cycles (estimated)
        cycles = []
        paths = []
        for i in range(len(weight_matrix)):
            for j in range(len(weight_matrix[0])):
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
            path1 = nx.shortest_path(graph, c[0], c[1], weight='weight')
            path2 = nx.shortest_path(graph, c[1], c[0], weight='weight')
            path = path1 + path2[1:]
            paths.append(path)
        return paths


class LibraryDesign(Design):

    def __init__(self, span_cost=None):
        super().__init__(span_cost)
        self.shared_alignments = []
        self._edges = []

    def _blast(self):
        super()._blast()
        self._share_query_blast()

    # @staticmethod
    # def _get_repeats_from_results(results):
    #     repeats = []
    #     for r in results:
    #         qk = r['query']['origin_key']
    #         sk = r['subject']['origin_key']
    #         if qk == sk:
    #             repeats.append((qk, r['query']['start'], r['query']['end']))
    #     return repeats

    def _get_repeats(self, alignments: List[Alignment]):
        for align in alignments:
            qk = align.query_key
            sk = align.subject_key
            if qk == sk:
                print("REPEAT")


    def _share_query_blast(self):
        self.logger.info("Finding shared fragments among queries")

        blast = self.blast_factory(self.QUERIES, self.QUERIES)
        blast.quick_blastn()
        results = blast.get_perfect()

        self.logger.info("Found {} shared alignments between the queries".format(len(results)))
        self.shared_alignments = results

        filtered_results = []
        self.container_factory.seqdb.update(blast.seq_db.records)
        self.container_factory.load_blast_json(results, Constants.SHARED_FRAGMENT)

        # TODO: expand the normal fragments with the shared fragments
        for container in self.container_list():
            # expand the share fragments using their own endpoints
            container.expand_pcr_products(container.get_groups_by_types(Constants.SHARED_FRAGMENT),
                                          Constants.SHARED_FRAGMENT)

            # expand the existing fragments with endpoints from the share alignments
            container.expand_pcr_products(container.get_groups_by_types(
                [Constants.FRAGMENT,
                Constants.PCR_PRODUCT,
                Constants.SHARED_FRAGMENT]
            ), Constants.PCR_PRODUCT)

            # grab the pcr products and expand primer pairs (again)
            templates = container.get_groups_by_types(
                Constants.PCR_PRODUCT
            )
            container.expand_primer_pairs(templates)


        repeats = []
        for container in self.container_list():
            # get all shared fragments
            alignments = container.get_alignments_by_types(Constants.SHARED_FRAGMENT)

            # add to list of possible repeats
            repeats += self._get_repeats(alignments)
        x = 1



        # # filter self alignmenbts
        #
        # for r in results:
        #     qk = r['query']['origin_key']
        #     sk = r['subject']['origin_key']
        #
        #     if qk != sk:
        #         n1 = (sk, r['subject']['start'], r['subject']['end'])
        #         n2 = (qk, r['query']['start'], r['query']['end'])
        #         if n1 not in repeats and n2 not in repeats:
        #             yield n1, n2

        # self.container_factory.load_blast_json(results, Constants.SHARED_FRAGMENT)
        # # TODO: need to expand certain points from these results...
        # # TODO: method that expands a list of points
        # # self.container_factory.load_blast_json(results, Constants.PCR_PRODUCT)
        #
        # edges = set()
        # for result in results:
        #     a = result['query']['start']
        #     b = result['query']['end']
        #     k = result['query']['origin_key']
        #     edges.add((a, b, k))
        # self._edges = edges
        # self.logger.info("{} possible ways to share fragments between {} goal plasmids.".format(2**len(edges), len(self.container_factory.alignments)))



    def optimize_library(self):
        combinations = range(2**len(self._edges))
        for e in self.logger.tqdm(combinations, "INFO", desc="optimizing for shared materials"):
            self.optimize()


        # # combine the sequence databases (for now)
        # seqdb = {}
        # seqdb.update(blast.seq_db.records)
        # seqdb.update(primer_blast.seq_db.records)
        #
        # self.container_factory = AlignmentContainerFactory(seqdb)
        # self.container_factory.load_blast_json(results, Constants.PCR_PRODUCT)
        # self.container_factory.load_blast_json(primer_results, Constants.PRIMER)