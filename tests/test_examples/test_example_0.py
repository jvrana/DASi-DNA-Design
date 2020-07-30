import json
import random
from os.path import abspath
from os.path import dirname
from os.path import join

import numpy as np
import pytest
from Bio import SeqIO

from dasi import Design
from dasi import LibraryDesign

here = abspath(dirname(__file__))
fixtures = join(here, "fixtures")


def mark(name, values, f=None, prefix="", suffix=""):
    if f is None:

        def f(x):
            return prefix + str(x) + suffix

    ids = [f(v) for v in values]
    return pytest.mark.parametrize(name, values, ids=ids)


@mark(
    "args",
    [
        # (3, 4, Design),
        (5, 6, Design),
        # (6, 7, Design),
        # (1, 2, Design),
        # (2, 3, Design),
        # (1, 6, Design),
        (None, None, LibraryDesign),
    ],
    prefix="index=",
)
def test_example_0(args):
    i0, i1, DesignCls = args
    random.seed(0)
    np.random.seed(0)

    def open_gb(path):
        with open(join(fixtures, path), "r") as f:
            return list(SeqIO.parse(f, format="genbank"))

    fragments = open_gb("fragments_0.gb")
    primers = open_gb("primers_0.gb")
    plasmids = open_gb("plasmids_0.gb")
    goals = open_gb("goals_0.gb")

    design = Design()
    design.add_fragments(fragments)
    design.add_primers(primers)
    design.add_templates(plasmids)
    s = slice(i0, i1)
    design.add_queries(goals[s])
    design.run(n_jobs=4)

    # assert successful runs
    print(design.status)
    for v in design.status.values():
        assert v["success"] is True

    # output JSON
    out = design.out()

    print(design.to_df()[1])
