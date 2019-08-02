"""Primer and synthesis design"""

from shoestring import AlignmentContainer, Constants, AssemblyGraphBuilder
from shoestring.utils import perfect_subject
import networkx as nx
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
from shoestring import BioBlastFactory
from shoestring.log import logger
from typing import List
from Bio.SeqRecord import SeqRecord


class Design(object):

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"

    def __init__(self):
        self.factory = BioBlastFactory()
        self.logger = logger(self)
        self.G = None

    def add(self, primers, templates, queries):
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

    def compile(self):
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

        container = AlignmentContainer(seqdb)

        # load the results (potential PCR Products)
        container.load_blast_json(results, Constants.PCR_PRODUCT)

        # load the primer results
        container.load_blast_json(primer_results, Constants.PRIMER)

        container.expand()

        # group by query_regions
        groups = container.alignment_groups

        self.logger.info("Number of types: {}".format(len(container.groups_by_type)))
        self.logger.info("Number of groups: {}".format(len(groups)))

        # build assembly graph
        graph_builder = AssemblyGraphBuilder(container)
        G = graph_builder.build_assembly_graph()

        self.logger.info("=== Assembly Graph ===")
        self.logger.info(nx.info(G))
        assert G.number_of_edges()
        self.G = G