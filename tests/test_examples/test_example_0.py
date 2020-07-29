import json
import random
from os.path import abspath
from os.path import dirname
from os.path import join

import numpy as np
import pytest
from Bio import SeqIO

from dasi import Design


here = abspath(dirname(__file__))
fixtures = join(here, "fixtures")


def mark(name, values, f=None, prefix="", suffix=""):
    if f is None:

        def f(x):
            return prefix + str(x) + suffix

    ids = [f(v) for v in values]
    return pytest.mark.parametrize(name, values, ids=ids)


@mark(
    "i",
    [
        (3, 4),
        # (5, 6),
        # (6, 7),
        # (1, 2),
        # (2, 3),
    ],
    prefix="index=",
)
def test_example_0(i):

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
    design.add_queries(goals[i[0] : i[1]])
    design.run()
    out = design.out()

    # print(json.dumps(out, indent=2))
    print(design.to_df()[1])
