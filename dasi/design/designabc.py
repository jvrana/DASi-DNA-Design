r"""
Design (:mod:`dasi.design`)
=============================

.. currentmodule:: dasi.design

This module provide DNA assembly functionality for DASi.
"""
import functools
import operator
from abc import ABC
from abc import abstractmethod
from multiprocessing import Pool
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import Tuple
from typing import TypeVar

import networkx as nx
import pandas as pd
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from pyblast import BioBlastFactory
from pyblast.utils import is_circular

from .design_result import DesignResult
from .graph_builder import AssemblyGraphBuilder
from .graph_builder import AssemblyGraphPostProcessor
from .optimize import optimize_graph
from dasi.__version__ import __title__
from dasi.__version__ import __version__
from dasi.constants import Constants
from dasi.cost import cached_span_cost
from dasi.cost import SpanCost
from dasi.design.graph_builder import AssemblyNode
from dasi.design.output import dasi_design_to_output_json
from dasi.design.output import validate_output
from dasi.design.report import Report
from dasi.exceptions import DasiDesignException
from dasi.log import logger
from dasi.models import Alignment
from dasi.models import AlignmentContainer
from dasi.models import AlignmentContainerFactory
from dasi.schemas import Schemas
from dasi.utils import chunkify
from dasi.utils import log_metadata
from dasi.utils import perfect_subject
from dasi.utils.sequence import generate_fake_designs

DesignType = TypeVar("Design")


BLAST_PENALTY_CONFIG = {
    "gapopen": 3,
    "gapextend": 3,
    "reward": 1,
    "penalty": -5,
    "ungapped": None,
}


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


class DesignABC(ABC):
    """Design class that returns optimal assemblies from a set of materials."""

    PRIMERS = "primers"
    TEMPLATES = "templates"
    QUERIES = "queries"
    FRAGMENTS = "fragments"
    DEFAULT_N_ASSEMBLIES = 3
    DEFAULT_JOB_SIZE = 1
    DEFAULT_N_JOBS = 1
    ALGORITHM = Constants.ALGORITHM_DEFAULT

    def __init__(self, span_cost=None, seqdb=None, name: str = ""):
        """

        :param span_cost:
        :type span_cost: SpanCost
        :param seqdb:
        :type seqdb: dict
        """
        self.name = name
        self.blast_factory = BioBlastFactory(config=BLAST_PENALTY_CONFIG)
        self.logger = logger(self)

        # graph by query_key
        if seqdb is None:
            seqdb = {}
        self._seqdb = seqdb  #: Sequence dict registry
        if span_cost is None:
            span_cost = cached_span_cost()
        self.span_cost = span_cost  #: span cost df
        self.graphs = {}  #: Dict[str, nx.DiGraph]
        self._results = {}  #: Dict[str, DesignResult]
        self.template_results = []
        self.fragment_results = []
        self.primer_results = []
        self.container_factory = AlignmentContainerFactory(self.seqdb)
        self._times = {}
        self._method_trace = {}

    def _get_design_status(self, qk):
        status = {"compiled": False, "run": False, "success": False, "assemblies": []}

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
            for a in self.results[qk].assemblies:
                status["success"] = True
                summ_df = a.to_df()
                material = sum(list(summ_df["material"]))
                eff = functools.reduce(operator.mul, summ_df["efficiency"])

                comp = 0
                for c in summ_df["complexity"]:
                    if not c:
                        continue
                    elif c > comp:
                        comp = c

                status["assemblies"].append(
                    {
                        "cost": {
                            "material cost": round(material, 2),
                            "assembly efficiency": round(eff, 2),
                            "max synthesis complexity": round(comp, 2),
                        }
                    }
                )
        return status

    @property
    def metadata(self):
        return {
            "program": __title__,
            "version": __version__,
            "execution_trace": self._method_trace,
        }

    @property
    def status(self):
        """Return the statuses for each subdesign."""
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
    @functools.wraps(generate_fake_designs)
    def fake(
        cls,
        n_designs: int,
        span_cost: SpanCost = None,
        circular: bool = True,
        n_linear_seqs: int = 50,
        n_cyclic_seqs: int = 50,
        n_primers: int = 50,
        n_primers_from_templates: int = 50,
        shared_length: int = 0,
        cyclic_size_int: Tuple[int, int] = (3000, 10000),
        linear_size_int: Tuple[int, int] = (100, 4000),
        primer_size_int: Tuple[int, int] = (15, 60),
        plasmid_size_interval: Tuple[int, int] = (5000, 10000),
        chunk_size_interval: Tuple[int, int] = (100, 3000),
        random_chunk_prob_int: Tuple[float, float] = (0, 0.5),
        random_chunk_size_int: Tuple[int, int] = (100, 1000),
        return_with_library: bool = False,
        **kwargs,
    ):
        library = generate_fake_designs(
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
            design_sequence_similarity_length=shared_length,
        )
        designs = library["design"]
        plasmids = library["cyclic"]
        fragments = library["linear"]
        primers = library["short"]

        design = cls(span_cost=span_cost)
        design.add_materials(
            primers=primers, fragments=fragments, templates=plasmids, queries=designs
        )
        if return_with_library:
            return design, library
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

    @property
    def container_list(self) -> List[AlignmentContainer]:
        """List of alignment containers in this design."""
        return list(self.container_factory.containers().values())

    @property
    def query_keys(self) -> List[str]:
        """List of query keys in this design."""
        return list(self.container_factory.containers())

    def assemble_graphs(self, query_keys: Optional[List[str]] = None):
        """Assemble all assembly graphs for all queries in this design."""
        for query_key, container in self.logger.tqdm(
            self._filter_containers(query_keys).items(),
            "INFO",
            desc="assembling graphs",
        ):
            self.graphs[query_key], _ = assemble_graph(container, self.span_cost)

    def post_process_graphs(self, post_processor_kwargs: Optional[dict] = None):
        if post_processor_kwargs is None:
            post_processor_kwargs = {}
        for qk, graph in self.graphs.items():
            container = self.container_factory.containers()[qk]
            query = self.seqdb[qk]
            processor = AssemblyGraphPostProcessor(
                graph,
                query,
                self.span_cost,
                self.seqdb,
                container,
                **post_processor_kwargs,
            )
            processor.logger = self.logger(processor)
            processor()

    def _filter_containers(self, query_keys):
        if query_keys is None:
            return self.container_factory.containers()
        return {
            k: v
            for k, v in self.container_factory.containers().items()
            if k in query_keys
        }

    def uncompile(self):
        self.container_factory.reset()
        self._results = {}

    @abstractmethod
    def precompile(self):
        pass

    @abstractmethod
    def postcompile(self, post_process_kwargs):
        pass

    @log_metadata("compile", additional_metadata={"algorithm": ALGORITHM})
    def compile(self, post_process_kwargs: Dict = None):
        """Compile the materials list into assembly graphs."""
        self.precompile()

        self.logger.info("Assembling graphs")
        self.assemble_graphs()
        self.postcompile(post_process_kwargs)

    def _run(self, n_paths: int, post_processing_kwargs: Dict):
        """Run the design. Runs `compile` and `optimize`, returning results
        that can be accessed by `design.results` or by `design.out()`

        :param n_paths: max number of assemblies per design to design
        :return: results
        """
        self.compile(post_processing_kwargs)
        return self.optimize(n_paths=n_paths)

    def run(
        self,
        n_paths: int = None,
        post_processing_kwargs: Dict = None,
        n_jobs: int = None,
        job_size: int = None,
    ):
        if n_paths is None:
            n_paths = self.DEFAULT_N_ASSEMBLIES
        if n_jobs is None:
            n_jobs = self.DEFAULT_N_JOBS
        if job_size is None:
            job_size = self.DEFAULT_JOB_SIZE
        if n_jobs > 1:
            self._run_with_pool(
                n_jobs=n_jobs,
                job_size=job_size,
                post_processing_kwargs=post_processing_kwargs,
            )
        else:
            self._run(n_paths=n_paths, post_processing_kwargs=post_processing_kwargs)
        self._freeze_graphs()

    # TODO: better multithreaded loggers
    @staticmethod
    def _pooled_run(args: List[Tuple[DesignType, str]]) -> Dict[str, DesignResult]:
        design, qk, n_paths, post_processing_kwargs, thread = args
        design.logger.name = "THREAD {}: ".format(thread) + design.logger.name

        design.assemble_graphs(query_keys=qk)
        design.postcompile(post_processing_kwargs)
        design.optimize(n_paths=n_paths)
        return design.graphs, design.results

    def _run_with_pool(
        self,
        n_jobs: int = None,
        job_size: int = None,
        n_paths: int = None,
        post_processing_kwargs: dict = None,
    ):
        if n_paths is None:
            n_paths = self.DEFAULT_N_ASSEMBLIES
        if n_jobs is None:
            n_jobs = self.DEFAULT_N_JOBS
        if job_size is None:
            job_size = self.DEFAULT_JOB_SIZE
        self.precompile()
        query_keys = self.query_keys
        if not query_keys:
            raise DasiDesignException("There are no design sequences.")

        # divide query_keys
        chunks = list(chunkify(query_keys, job_size))
        print(chunks)
        args = list(
            zip(
                [self] * len(chunks),
                chunks,
                [n_paths] * len(chunks),
                [post_processing_kwargs] * len(chunks),
                range(len(chunks)),
            )
        )
        with Pool(processes=min(len(chunks), n_jobs)) as pool:
            graphs_and_results = pool.map(self._pooled_run, args)

        # update graphs and results
        combined_results = {}
        combined_graphs = {}
        for g, r in graphs_and_results:
            combined_results.update(r)
            combined_graphs.update(g)
        self.graphs = combined_graphs
        self._results = combined_results
        return self.results

    @property
    def is_compiled(self):
        """Return whether the design has been compiled.

        :return: True if compiled. False if otherwise.
        """
        if self.container_list:
            return True
        return False

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

    @log_metadata("optimize", additional_metadata={"algorithm": ALGORITHM})
    def optimize(self, n_paths=DEFAULT_N_ASSEMBLIES) -> Dict[str, DesignResult]:

        if not self.container_list:
            raise DasiDesignException(
                "Design must be compiled before running optimization.'"
            )
        self._results = self._optimize(n_paths)
        return self._results

    def _freeze_graphs(self):
        for graph in self.graphs.values():
            nx.freeze(graph)

    def _optimize(self, n_paths) -> Dict[str, DesignResult]:
        """Finds the optimal paths for each query in the design."""
        results_dict = {}
        for query_key, graph, query_length, cyclic, result in self.logger.tqdm(
            self._collect_optimize_args(self.graphs),
            "INFO",
            desc="optimizing graphs (n_graphs={})".format(len(self.graphs)),
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

    def report(self):
        return Report(self)

    @property
    def last_run_start(self):
        x = self._method_run_times.get("optimize", None)
        if not x:
            return None
        else:
            return x[0]

    @property
    def last_run_end(self):
        x = self._method_run_times.get("optimize", None)
        if not x:
            return None
        else:
            return x[1]

    @property
    def last_compile_start(self):
        x = self._method_run_times.get("compile", None)
        if not x:
            return None
        else:
            return x[0]

    @property
    def last_compile_end(self):
        x = self._method_run_times.get("compile", None)
        if not x:
            return None
        else:
            return x[1]

    def out(
        self,
        fmt: str = "json",
        elim_extra_reactions: bool = False,
        query_keys: Optional[List[str]] = None,
    ):
        """Return the results of the design as a validates output JSON.

        The output JSON is follows the following schema, see
        :param fmt:
        :return:
        """
        if fmt.lower() == "json":
            output = dasi_design_to_output_json(
                self, elim_extra_reactions=elim_extra_reactions, query_keys=query_keys
            )
            validate_output(output)
            return output
        else:
            raise ValueError("Format '{}' not recognized".format(fmt))

    def out_jsonschema(self):
        """Return the output JSONschema."""
        return Schemas.output_schema
