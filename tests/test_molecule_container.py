from shoestring import AlignmentContainer, Constants
from shoestring.utils import perfect_subject
import networkx as nx
from os.path import join


def test_produce_assembly_graph(blast_factory, here):
    blast = blast_factory("templates", "queries")
    primer_blast = blast_factory("primers", "queries")

    blast.quick_blastn()
    primer_blast.quick_blastn_short()

    results = blast.get_perfect()
    primer_results = primer_blast.get_perfect()

    primer_results = [p for p in primer_results if perfect_subject(p["subject"])]
    print("Number of perfect primers: {}".format(len(primer_results)))
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

    # group by query_regions
    groups = container.alignment_groups

    print("Number of types: {}".format(len(container.groups_by_type)))
    print("Number of groups: {}".format(len(groups)))

    # build assembly graph
    G = container.build_assembly_graph()
    # print()
    print("=== Assembly Graph ===")
    print(nx.info(G))
    assert G.number_of_edges()

    # print_graph
    nx.write_gexf(G, join(here, "out", "graph.gexf"))
