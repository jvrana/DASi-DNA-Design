from multiprocessing import Pool
from typing import List
from typing import Tuple

import networkx as nx

from .optimize import optimize_graph
from dasi.cost import SpanCost
from dasi.design.graph_builder import AssemblyGraphBuilder
from dasi.models import AlignmentContainer
from dasi.models import AlignmentContainerFactory


def _multiprocessing_optimize_graph(
    args: Tuple[nx.DiGraph, int, bool, int]
) -> List[List[tuple]]:
    paths, _ = optimize_graph(args[0], args[1], args[2], args[3])
    return paths


def multiprocessing_optimize_graph(
    graphs: List[nx.DiGraph],
    query_lengths: List[int],
    cyclics: List[bool],
    n_paths: int,
    n_jobs: int,
):
    """Optimize graphs using multiprocessing."""
    args = [(g, q, c, n_paths) for g, q, c in zip(graphs, query_lengths, cyclics)]

    with Pool(processes=min(n_jobs, len(graphs))) as pool:  # start 4 worker processes
        paths = pool.map(_multiprocessing_optimize_graph, args)
    return paths


def assemble_graph(
    container: AlignmentContainer, span_cost: SpanCost
) -> Tuple[nx.DiGraph, AlignmentContainer]:
    """Build an assembly graph for a specified query."""
    container.expand(expand_overlaps=True, expand_primers=True)
    container.clean_alignments()
    container.groups()
    container.freeze()
    graph_builder = AssemblyGraphBuilder(container, span_cost=span_cost)
    graph = graph_builder.build_assembly_graph()
    return graph, container


def _multiprocessing_assemble_graph(
    arg: Tuple[AlignmentContainer, SpanCost]
) -> Tuple[nx.DiGraph, AlignmentContainer]:
    return assemble_graph(arg[0], arg[1])


def multiprocessing_assemble_graph(
    container_factory: AlignmentContainerFactory, span_cost: SpanCost, n_jobs: int
) -> List[nx.DiGraph]:
    """Assemble graphs using multiprocessing."""
    query_keys = container_factory.alignments.keys()
    containers = [container_factory.containers()[k] for k in query_keys]

    args = [(container, span_cost) for container in containers]
    with Pool(
        processes=min(n_jobs, len(containers))
    ) as pool:  # start 4 worker processes
        graphs, expanded_containers = zip(
            *pool.map(_multiprocessing_assemble_graph, args)
        )

    # update container_factory alignments

    new_containers = dict(zip(query_keys, expanded_containers))
    # new_alignments = {key: c.alignments for key, c in new_containers.items()}
    # container_factory.set_alignments
    for key, container in container_factory.containers().items():
        container_factory._alignments[key] = new_containers[key].alignments
        container_factory._containers[key] = new_containers[key]
        container.seqdb = container_factory.seqdb
    return graphs
