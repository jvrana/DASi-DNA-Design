import random
from multiprocessing import Pool

import numpy as np
import pytest

from dasi.design import Design
from dasi.design import LibraryDesign


@pytest.mark.parametrize("i", [1, 2, 3])
def test_compile_with_filtered_keys(i):
    design = Design.fake(3)
    design._blast()
    query_keys = design.query_keys
    design.assemble_graphs(query_keys=query_keys[:i])
    design.post_process_graphs(**{})
    design.optimize()
    assert len(design.results) == i


@pytest.mark.parametrize("n", [1, 16])
@pytest.mark.parametrize("nseqs", [3, 1])
@pytest.mark.parametrize("s", [1, 2])
def test_pooled_run(n, nseqs, s):
    design = LibraryDesign.fake(nseqs)
    design._run_with_pool(n_jobs=n, job_size=s)
    assert len(design.graphs) == nseqs
    assert len(design.results) == nseqs


def test_compare_pooled_run():
    """Results of pooled and non-pooled ought to be the same."""

    def pooled_run():
        random.seed(0)
        np.random.seed(0)

        design = LibraryDesign.fake(10)
        design._run_with_pool(16, 1)
        return design

    def run():
        random.seed(0)
        np.random.seed(0)

        design = LibraryDesign.fake(10)
        design.run()
        return design

    design1 = run()
    design2 = pooled_run()
    df1 = design1.to_df()[1]
    df2 = design2.to_df()[1]

    keys = [
        "query_start",
        "query_end",
        "cost",
        "material",
        "span",
        "type",
        "internal_or_external",
        "efficiency",
        "complexity",
        "notes",
        "DESIGN_ID",
    ]
    assert df1[keys].equals(df2[keys])
