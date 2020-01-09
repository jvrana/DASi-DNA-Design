r"""
Command Line (:mod:`dasi.command_line`)
=============================

.. currentmodule:: dasi.command_line

This module provide command line interface for DASi.

.. code-block::

    dasi run
"""
import os

import fire
from Bio import BiopythonParserWarning
from Bio import SeqIO
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import __version__
from dasi import Design
from dasi.cost import SpanCost
from dasi.log import logger


class DasiCLI:
    """DASi command line interface."""

    def __init__(
        self,
        directory=os.getcwd(),
        primers="primers.fasta",
        templates="templates/*.gb",
        fragments="fragments/*.gb",
        goals="goals/*.gb",
        verbose="v",
    ):
        """Initialize a new design."""
        self._directory = directory
        self._primers = os.path.join(self._directory, primers)
        self._templates = os.path.join(self._directory, templates)
        self._fragments = os.path.join(self._directory, fragments)
        self._goals = os.path.join(self._directory, goals)
        self._do_save = True
        self._logger = logger(self)
        if verbose == "v":
            self._logger.set_level("INFO")
        elif verbose == "vv":
            logger.set_level("INFO")
        elif verbose == "vvv":
            logger.set_level("DEBUG")
        elif verbose is None:
            pass
        else:
            raise ValueError(
                "Verbose level '{}' not recognized. " "Select from 'v', 'vv', or 'vvv'"
            )

    def run(self, n_jobs: int = 10):
        """Run a design job.

        :param n_jobs: number of parrallel jobs to run. (default: 10)
        :return:
        """
        import warnings

        warnings.simplefilter(action="ignore", category=RuntimeWarning)
        warnings.simplefilter(action="ignore", category=BiopythonParserWarning)

        self._logger.info("Loading sequence files")
        primers = make_linear(load_fasta_glob(self._primers))
        templates = make_circular(load_genbank_glob(self._templates))
        fragments = make_linear(load_genbank_glob(self._fragments))
        goals = make_circular(load_genbank_glob(self._goals))
        design = Design()
        design.n_jobs = n_jobs
        design.add_materials(
            primers=primers, templates=templates, fragments=fragments, queries=goals
        )

        self._logger.info("Getting span cost model")
        span_cost = self._get_span_cost()
        design.span_cost = span_cost

        self._logger.info("Compiling possible molecular assemblies")
        design.compile()

        self._logger.info("Optimizing molecular assemblies")
        design.optimize()

        self._logger.info("Designing assembly primers and fragments")
        df, adf, design_json = design.to_df()
        adf.to_csv("summary.csv")
        df.to_csv("sequence_design.csv")

        records = []
        for result in design.results.values():
            if result.assemblies:
                a = result.assemblies[0]
                for i, role, m in a.molecules:
                    records.append(m.sequence)

        SeqIO.write(records, os.path.join(self._directory, "sequences.gb"), "genbank")

    def version(self):
        """Print the package version."""
        print(__version__)

    def _get_span_cost(self):
        """Saves the span cost as bytes; reloads when called."""
        path = os.path.join(self._directory, "span_cost.b")
        if self._do_save and os.path.isfile(path):
            with logger.timeit("INFO", "loading bytes"):
                print("Loading file: {}".format(path))
                span_cost = SpanCost.load(path)
        else:
            span_cost = SpanCost.open()
            if self._do_save:
                with logger.timeit("INFO", "saving bytes"):
                    print("Saving file: {}".format(path))
                    span_cost.dump(path)
        return span_cost


def main():
    fire.Fire(DasiCLI)
