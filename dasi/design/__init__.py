r"""
Design (:mod:`dasi.design`)
=============================

.. currentmodule:: dasi.design

This module provide DNA assembly functionality for DASi.
"""
from __future__ import annotations

import bisect
from collections.abc import Iterable
from typing import Dict
from typing import Generator
from typing import List
from typing import Tuple

import networkx as nx
import pandas as pd
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from pyblast import BioBlastFactory
from pyblast.utils import is_circular

from .design_algorithms import assemble_graph
from .design_algorithms import multiprocessing_assemble_graph
from .design_algorithms import multiprocessing_optimize_graph
from .design_algorithms import optimize_graph
from .graph_builder import AssemblyGraphBuilder
from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.design.graph_builder import AssemblyNode
from dasi.exceptions import DasiInvalidMolecularAssembly
from dasi.exceptions import DASiWarning
from dasi.log import logger
from dasi.models import Alignment
from dasi.models import AlignmentContainer
from dasi.models import AlignmentContainerFactory
from dasi.models import Assembly
from dasi.utils import perfect_subject
from dasi.utils import prep_df

BLAST_PENALTY_CONFIG = {"gapopen": 3, "gapextend": 3, "reward": 1, "penalty": -5}


class DesignResult(Iterable):
    """DesignResult container.

    Maintains a list of top assemblies.
    """

    def __init__(
        self, container: AlignmentContainer, graph: nx.DiGraph, query_key: str
    ):
        self.container = container
        self.graph = graph
        self.query_key = query_key
        self.query = self.container.seqdb[query_key]
        self._assemblies = []
        self._keys = []

    @property
    def seqdb(self):
        return self.container.seqdb

    @property
    def assemblies(self) -> Tuple[Assembly, ...]:
        """Return a tuple of all assemblies.

        :return: tuple of all assemblies.
        """
        return tuple(self._assemblies)

    def _add_assembly_from_path(self, path: List[AssemblyNode]):
        return Assembly(
            path,
            self.container,
            self.graph,
            self.query_key,
            self.query,
            seqdb=self.seqdb,
            do_raise=False,
        )

    def add_assembly(
        self,
        path: List[AssemblyNode],
        ignore_invalid: bool = False,
        allow_invalid: bool = False,
    ):
        """Add an assembly from a list of nodes.

        :param path: list of nodes
        :return: None
        """

        assembly = self._add_assembly_from_path(path)

        try:
            assembly.post_validate()
        except DasiInvalidMolecularAssembly as e:
            if ignore_invalid:
                return
            elif not allow_invalid:
                raise e

        cost = assembly.cost()
        n_nodes = len(assembly._nodes)
        k = (cost, n_nodes)
        i = bisect.bisect_left(self._keys, k)
        self._assemblies.insert(i, assembly)
        self._keys.insert(i, k)
        return assembly

    def add_assemblies(
        self,
        paths: List[List[AssemblyNode]],
        ignore_invalid: bool = False,
        allow_invalid: bool = False,
    ):
        """Adds a list of assemblies.

        :param paths: list of list of paths
        :return: None
        """
        for path in paths:
            self.add_assembly(
                path, ignore_invalid=ignore_invalid, allow_invalid=allow_invalid
            )

    #
    #
    #     for a in self.assemblies:
    #         for n1, n2, edata in a.edges():
    #             result = edata["sequence_result"]
    #             mol_type = edata["type_def"]
    #             if mol_type.int_or_ext == "internal":
    #                 pass
    #                 # add a new reaction
    #                 # add template
    #                 # add primers
    #             elif mol_type.use_direct:
    #                 pass
    #                 # add fragment
    #             elif mol_type.synthesize:
    #                 pass
    #                 # add synthesize
    #             else:
    #                 pass
    #                 # raise Exception

    def __iter__(self) -> Generator[Assembly]:
        """Yield assemblies.

        :yield: assembly
        """
        for assembly in self.assemblies:
            yield assembly

    def __getitem__(self, item: str) -> Assembly:
        return list(self)[item]

    def __str__(self):
        return "<{cls} query={qname} {qk} nassemblies={n}>".format(
            cls=self.__class__.__name__,
            qname=self.query.name,
            qk=self.query_key,
            n=len(self.assemblies),
        )


class FakePool:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return FakePool()

    def __exit__(self, a, b, c):
        pass

    @staticmethod
    def map(func, args):
        return [func(arg) for arg in args]


class Design:
    """Design class that returns optimal assemblies from a set of materials."""

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"
    FRAGMENTS = "fragments"
    DEFAULT_N_JOBS = 1

    def __init__(self, span_cost=None, seqdb=None, n_jobs=None):
        """

        :param span_cost:
        :type span_cost: SpanCost
        :param seqdb:
        :type seqdb: dict
        :param n_jobs:
        :type n_jobs: int
        """
        self.blast_factory = BioBlastFactory()
        self.logger = logger(self)

        # graph by query_key
        if seqdb is None:
            seqdb = {}
        self._seqdb = seqdb  #: Sequence dict registry
        self.span_cost = span_cost  #: span cost df
        self.graphs = {}  #: Dict[str, nx.DiGraph]
        self._results = {}  #: Dict[str, DesignResult]
        self.template_results = []
        self.fragment_results = []
        self.primer_results = []
        self.container_factory = AlignmentContainerFactory(self.seqdb)
        self.n_jobs = n_jobs or self.DEFAULT_N_JOBS  #: number of multiprocessing jobs

    @property
    def seqdb(self) -> Dict[str, SeqRecord]:
        return self._seqdb

    @property
    def results(self):
        return dict(self._results)

    def add_materials(
        self,
        primers: List[SeqRecord],
        templates: List[SeqRecord],
        queries: List[SeqRecord],
        fragments: List[SeqRecord] = None,
    ):
        if fragments is None:
            fragments = []
        self.add_primers(primers)
        fragments = self.filter_linear_records(fragments)
        self.add_templates(templates + fragments)
        self.add_queries(queries)
        self.add_fragments(fragments)

    def add_primers(self, primers: List[SeqRecord]):
        """Add primer sequences to materials list."""
        self.blast_factory.add_records(primers, self.PRIMERS)
        self.logger.info("Added {} primers".format(len(primers)))

    def add_templates(self, templates: List[SeqRecord]):
        """Add template sequences to materials list."""
        self.blast_factory.add_records(templates, self.TEMPLATES)
        self.logger.info("Added {} templates".format(len(templates)))

    def add_queries(self, queries: List[SeqRecord]):
        """Add goal/query sequences to materials list."""
        self.blast_factory.add_records(queries, self.QUERIES)
        self.logger.info("Added {} queries".format(len(queries)))

    def add_fragments(self, fragments: List[SeqRecord]):
        """Add fragment sequences to materials list."""
        self.blast_factory.add_records(fragments, self.FRAGMENTS)
        self.logger.info("Added {} queries".format(len(fragments)))
        # self.blast_factory.add_records(fragments, self.TEMPLATES)

    @classmethod
    def filter_linear_records(cls, records: List[SeqRecord]) -> List[SeqRecord]:
        """Return only linear records."""
        return [r for r in records if not is_circular(r)]

    @classmethod
    def filter_perfect_subject(cls, results: dict) -> List[dict]:
        """return only results whose subject is 100% aligned to query."""
        return [r for r in results if perfect_subject(r["subject"])]

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
    def containers(self) -> Dict[str, AlignmentContainer]:
        """Iterable of alignment containers in this design."""
        return self.container_factory.containers()

    def container_list(self) -> List[AlignmentContainer]:
        """List of alignment containers in this design."""
        return list(self.container_factory.containers().values())

    def query_keys(self) -> List[str]:
        """List of query keys in this design."""
        return list(self.container_factory.containers())

    def assemble_graphs(self, n_jobs=None):
        n_jobs = n_jobs or self.n_jobs
        if n_jobs > 1:
            with self.logger.timeit(
                "DEBUG",
                "assembling graphs (n_graphs={}, threads={})".format(
                    len(self.container_list()), n_jobs
                ),
            ):
                self._assemble_graphs_with_threads(n_jobs)
        else:
            self._assemble_graphs_without_threads()

    def _assemble_graphs_without_threads(self):
        """Assemble all assembly graphs for all queries in this design."""
        for query_key, container in self.logger.tqdm(
            self.container_factory.containers().items(),
            "INFO",
            desc="assembling graphs (threads=1)",
        ):
            self.graphs[query_key], _ = assemble_graph(container, self.span_cost)

    def _assemble_graphs_with_threads(self, n_jobs=None):
        query_keys, containers = zip(*self.container_factory.containers().items())

        graphs = multiprocessing_assemble_graph(
            self.container_factory, self.span_cost, n_jobs=n_jobs
        )

        # update graphs dict
        for qk, g, c in zip(query_keys, graphs, containers):
            self.graphs[qk] = g

    def compile(self, n_jobs=None):
        """Compile materials to assembly graph."""
        self._results = {}
        with self.logger.timeit("DEBUG", "running blast"):
            self._blast()
        self.assemble_graphs(n_jobs=n_jobs)

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
    def _find_iter_alignment(a: int, b: int, alignments: Iterable[Alignment]):
        for align in alignments:
            if a == align.query_region.a and b == align.query_region.b:
                yield align

    @staticmethod
    def path_to_edge_costs(
        path: List[AssemblyNode], graph: nx.DiGraph
    ) -> List[Tuple[AssemblyNode, AssemblyNode, dict]]:
        arr = []
        for n1, n2 in pairwise(path):
            edata = graph[n1][n2]
            arr.append((n1, n2, edata))
        return arr

    def optimize(self, n_paths=3, n_jobs=None) -> Dict[str, DesignResult]:
        n_jobs = n_jobs or self.n_jobs
        if n_jobs > 1:
            with self.logger.timeit(
                "DEBUG",
                "optimizing graphs (n_graphs={}, threads={})".format(
                    len(self.graphs), n_jobs
                ),
            ):
                self._results = self._optimize_with_threads(n_paths, n_jobs)
        else:
            self._results = self._optimize_without_threads(n_paths)
        return self._results

    def _optimize_without_threads(self, n_paths) -> Dict[str, DesignResult]:
        """Finds the optimal paths for each query in the design."""
        results_dict = {}
        for query_key, graph, query_length, cyclic, result in self.logger.tqdm(
            self._collect_optimize_args(self.graphs),
            "INFO",
            desc="optimizing graphs (n_graphs={}, threads=1)".format(len(self.graphs)),
        ):

            container = self.containers[query_key]
            query = container.seqdb[query_key]
            cyclic = is_circular(query)
            results_dict[query_key] = result
            paths, costs = optimize_graph(graph, len(query), cyclic, n_paths)
            if not paths:
                query_rec = self.blast_factory.db.records[query_key]
                self.logger.error(
                    "\n\tThere were no solutions found for design '{}' ({}).\n\t"
                    "This sequence may be better synthesized. Use a tool such as JBEI's"
                    " BOOST.".format(query_rec.name, query_key)
                )
            result.add_assemblies(paths, ignore_invalid=True)
        return results_dict

    def _optimize_with_threads(self, n_paths=5, n_jobs=10) -> Dict[str, DesignResult]:
        results_dict = {}
        query_keys, graphs, query_lengths, cyclics, results = zip(
            *list(self._collect_optimize_args(self.graphs))
        )

        list_of_paths = multiprocessing_optimize_graph(
            graphs=graphs,
            query_lengths=query_lengths,
            cyclics=cyclics,
            n_paths=n_paths,
            n_jobs=n_jobs,
        )
        for qk, paths, result in zip(query_keys, list_of_paths, results):
            result.add_assemblies(paths, ignore_invalid=True)
            results_dict[qk] = result
        return results_dict

    def _collect_optimize_args(
        self, graphs: Dict[str, nx.DiGraph]
    ) -> Tuple[str, nx.DiGraph, bool, dict]:
        for query_key, graph in self.logger.tqdm(
            graphs.items(), "INFO", desc="optimizing graphs"
        ):
            container = self.containers[query_key]
            query = container.seqdb[query_key]
            cyclic = is_circular(query)
            result = DesignResult(container=container, query_key=query_key, graph=graph)
            yield query_key, graph, len(query), cyclic, result

    # TODO: order keys
    def to_df(self, assembly_index=0):
        dfs = []
        adfs = []
        for i, (qk, result) in enumerate(self.results.items()):
            if result.assemblies:
                a = result.assemblies[assembly_index]
                df = a.to_reaction_df()
                df["DESIGN_ID"] = i
                df["DESIGN_KEY"] = qk
                df["ASSEMBLY_ID"] = assembly_index
                dfs.append(df)

                adf = a.to_df()
                adf["DESIGN_ID"] = i
                adf["DESIGN_KEY"] = qk
                adf["ASSEMBLY_ID"] = assembly_index
                adfs.append(adf)
            else:
                msg = "Query {} {} yielded no assembly".format(self.seqdb[qk].name, i)
                self.logger.error(msg)

        df = pd.concat(dfs)
        colnames = [
            "DESIGN_ID",
            "DESIGN_KEY",
            "ASSEMBLY_ID",
            "REACTION_ID",
            "NAME",
            "TYPE",
            "KEY",
            "ROLE",
            "REGION",
            "SEQUENCE",
            "LENGTH",
            "META",
        ]
        df.columns = colnames
        df.sort_values(
            by=["TYPE", "DESIGN_ID", "ASSEMBLY_ID", "REACTION_ID", "NAME", "ROLE"],
            inplace=True,
        )

        adf = pd.concat(adfs)
        return df, adf


class LibraryDesign(Design):
    """Design class for producing assemblies for libraries."""

    DEFAULT_N_JOBS = 10

    def __init__(self, span_cost=None, n_jobs=None):
        super().__init__(span_cost=span_cost, n_jobs=n_jobs)
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

    # TODO: why?
    @staticmethod
    def _get_iter_non_repeats(alignments: List[Alignment]) -> Tuple[str, int, int]:
        """Return repeat regions of alignments. These are alignments that align
        to themselves.

        :param alignments:
        :return:
        """
        for align in alignments:
            qk = align.query_key
            sk = align.subject_key
            if qk == sk:
                yield (qk, align.query_region.a, align.query_region.b)

    def _share_query_blast(self):
        """Find and use shared fragments across queries.

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
            non_repeats = list(self._get_iter_non_repeats(alignments))
            self.logger.info(
                "{} non repeats for {}".format(len(non_repeats), query_key)
            )
            # add to list of possible repeats
            # repeats += list(self._get_iter_non_repeats(alignments))
        self.repeats = repeats

    def compile_library(self, n_jobs=None):
        """Compile the materials list into assembly graphs."""
        n_jobs = n_jobs or self.DEFAULT_N_JOBS
        self.graphs = {}
        self._blast()
        self._share_query_blast()
        self.assemble_graphs(n_jobs=n_jobs)

    def optimize_library(self):
        """Optimize the assembly graph for library assembly."""
        raise NotImplementedError
