r"""
Design (:mod:`dasi.design`)
=============================

.. currentmodule:: dasi.design

This module provide DNA assembly functionality for DASi.
"""
import bisect
import functools
import operator
from typing import Dict
from typing import Generator
from typing import Iterable
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
from .graph_builder import AssemblyGraphPostProcessor
from .optimize import optimize_graph
from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.design.graph_builder import AssemblyNode
from dasi.exceptions import DasiInvalidMolecularAssembly
from dasi.log import logger
from dasi.models import Alignment
from dasi.models import AlignmentContainer
from dasi.models import AlignmentContainerFactory
from dasi.models import AlignmentGroup
from dasi.models import Assembly
from dasi.utils import perfect_subject
from dasi.utils.testing_utils import fake_designs

BLAST_PENALTY_CONFIG = {
    "gapopen": 3,
    "gapextend": 3,
    "reward": 1,
    "penalty": -5,
    "ungapped": None,
}


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

    def __iter__(self) -> Generator[Assembly, None, None]:
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
        self.blast_factory = BioBlastFactory(config=BLAST_PENALTY_CONFIG)
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

    def _get_design_status(self, qk):
        status = {
            "compiled": False,
            "run": False,
            "failed": False,
            "success": False,
            "cost": {},
        }

        record = self.seqdb[qk]
        status["record"] = {
            "name": record.name,
            "length": len(record.seq),
            "id": record.id,
            "is_circular": is_circular(record),
        }

        if self.graphs.get(qk, None) is not None:
            status["compiled"] = True

        if self.results.get(qk, None) is not None:
            status["run"] = True
            if self.results[qk].assemblies:
                status["success"] = True
                summ_df = self.results[qk].assemblies[0].to_df()
                material = sum(list(summ_df["material"]))
                eff = functools.reduce(operator.mul, summ_df["efficiency"])

                comp = 0
                for c in summ_df["complexity"]:
                    if not c:
                        continue
                    elif c > comp:
                        comp = c

                status["cost"] = {
                    "material": round(material, 2),
                    "efficiency": round(eff, 2),
                    "max_synth_complexity": round(comp, 2),
                }
            else:
                status["failed"] = True
        return status

    @property
    def status(self):
        statuses = {}
        for qk in self.container_factory.containers():
            statuses[qk] = self._get_design_status(qk)
        return statuses

    @property
    def seqdb(self) -> Dict[str, SeqRecord]:
        return self._seqdb

    @property
    def results(self):
        return dict(self._results)

    @classmethod
    @functools.wraps(fake_designs)
    def fake(
        cls,
        n_designs: int,
        span_cost: SpanCost = None,
        circular: bool = True,
        n_linear_seqs: int = 50,
        n_cyclic_seqs: int = 50,
        n_primers: int = 50,
        n_primers_from_templates: int = 50,
        cyclic_size_int: Tuple[int, int] = (3000, 10000),
        linear_size_int: Tuple[int, int] = (100, 4000),
        primer_size_int: Tuple[int, int] = (15, 60),
        plasmid_size_interval: Tuple[int, int] = (5000, 10000),
        chunk_size_interval: Tuple[int, int] = (100, 3000),
        random_chunk_prob_int: Tuple[float, float] = (0, 0.5),
        random_chunk_size_int: Tuple[int, int] = (100, 1000),
        **kwargs,
    ):
        library = fake_designs(
            n_designs=n_designs,
            circular=circular,
            n_linear_seqs=n_linear_seqs,
            n_cyclic_seqs=n_cyclic_seqs,
            n_primers=n_primers,
            n_primers_from_templates=n_primers_from_templates,
            cyclic_size_int=cyclic_size_int,
            linear_size_int=linear_size_int,
            primer_size_int=primer_size_int,
            plasmid_size_interval=plasmid_size_interval,
            chunk_size_interval=chunk_size_interval,
            random_chunk_prob_int=random_chunk_prob_int,
            random_chunk_size_int=random_chunk_size_int,
            **kwargs,
        )
        designs = library["design"]
        plasmids = library["cyclic"]
        fragments = library["linear"]
        primers = library["short"]

        design = cls(span_cost=span_cost)
        design.add_materials(
            primers=primers, fragments=fragments, templates=plasmids, queries=designs
        )
        return design

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
        template_blast.blastn()
        template_results = template_blast.get_perfect()
        self.template_results = template_results
        self.logger.info("Number of template matches: {}".format(len(template_results)))
        self.container_factory.seqdb.update(template_blast.seq_db.records)

        # align fragments
        if self.blast_factory.record_groups[self.FRAGMENTS]:
            fragment_blast = self.blast_factory(self.FRAGMENTS, self.QUERIES)
            fragment_blast.blastn()
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
            primer_blast.blastn_short()
            primer_results = self.filter_perfect_subject(primer_blast.get_perfect())
            self.container_factory.seqdb.update(primer_blast.seq_db.records)
            self.logger.info(
                "Number of perfect primers: {}".format(len(primer_results))
            )
        else:
            primer_results = []
        self.primer_results = primer_results

        # initialize container factory
        query_keys = [rec.id for rec in self.blast_factory.record_groups[self.QUERIES]]
        query_keys = [self.blast_factory.db.get_origin_key(k) for k in query_keys]
        self.container_factory.initialize_empty(query_keys)

        # load results
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

    def post_process_graphs(self):
        for qk, graph in self.graphs.items():
            query = self.seqdb[qk]
            processor = AssemblyGraphPostProcessor(graph, query, self.span_cost)
            processor()

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
        self.post_process_graphs()

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

        self._freeze_graphs()
        return self._results

    def _freeze_graphs(self):
        for graph in self.graphs.values():
            nx.freeze(graph)

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
    # TODO: group identical reactions (identical output sequence)
    def to_df(self, assembly_index: int = 0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        reaction_dfs = []
        summary_dfs = []
        design_json = {}
        for i, (qk, result) in enumerate(self.results.items()):
            if result.assemblies and result.assemblies[0]:
                assembly = result.assemblies[assembly_index]
                react_df = assembly.to_reaction_df()
                react_df["DESIGN_ID"] = i
                react_df["DESIGN_KEY"] = qk
                react_df["ASSEMBLY_ID"] = assembly_index
                reaction_dfs.append(react_df)

                summ_df = assembly.to_df()
                summ_df["DESIGN_ID"] = i
                summ_df["DESIGN_KEY"] = qk
                summ_df["ASSEMBLY_ID"] = assembly_index
                summary_dfs.append(summ_df)

                rec = self.seqdb[qk]
                design_json[qk] = {"NAME": rec.name, "FILE": None}
                if hasattr(rec, "from_file"):
                    design_json[qk]["FILE"] = rec.from_file
            else:
                msg = "Query {} {} yielded no assembly".format(self.seqdb[qk].name, i)
                self.logger.error(msg)

        colnames = [
            "DESIGN_ID",
            "DESIGN_KEY",
            "ASSEMBLY_ID",
            "REACTION_ID",
            "REACTION_NAME",
            "NAME",
            "TYPE",
            "KEY",
            "ROLE",
            "REGION",
            "SEQUENCE",
            "LENGTH",
            "META",
        ]
        if reaction_dfs:
            react_df = pd.concat(reaction_dfs)
        else:
            react_df = pd.DataFrame(columns=colnames)
        react_df.columns = colnames
        react_df.sort_values(
            by=["TYPE", "DESIGN_ID", "ASSEMBLY_ID", "REACTION_ID", "NAME", "ROLE"],
            inplace=True,
        )

        if summary_dfs:
            summ_df = pd.concat(summary_dfs)
        else:
            summ_df = pd.DataFrame()

        return react_df, summ_df, design_json
