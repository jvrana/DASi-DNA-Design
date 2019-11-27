from copy import deepcopy
from os.path import join

import networkx as nx
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.constants import Constants
from dasi.design import LibraryDesign
from dasi.design.graph_builder import AssemblyGraphPostProcessor
from dasi.models import AlignmentContainer
from dasi.models import AlignmentGroup
from dasi.utils import group_by
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


def test_library_design_draft(paths, here, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))
    queries = queries

    design = LibraryDesign(span_cost=span_cost)
    design.n_jobs = 1
    design.add_materials(primers=primers, templates=templates, queries=queries)
    design.compile_library()
    results = design.optimize_library()

    for qk, result in results.items():
        df = result.assemblies[0].to_df()
        print(design.seqdb[qk].name)
        print(df)
