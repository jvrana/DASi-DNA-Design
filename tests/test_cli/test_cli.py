import os
import shutil

import pytest

from dasi.cli import CLI

here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture
def outdir(paths):
    out = os.path.join(here, ".test_out")
    if not os.path.isdir(out):
        os.mkdir(out)
    templates = os.path.join(out, "templates")
    goals = os.path.join(out, "goals")
    primers = os.path.join(out, "primers.fasta")

    if os.path.isdir(templates):
        shutil.rmtree(templates)

    if os.path.isdir(goals):
        shutil.rmtree(goals)

    if os.path.isfile(primers):
        os.remove(primers)

    shutil.copytree(os.path.dirname(paths["registry"]), templates)
    shutil.copytree(os.path.dirname(paths["queries"]), goals)
    shutil.copy(paths["primers"], os.path.join(out, primers))
    yield out
    # shutil.rmtree(os.path.join(out, 'templates'))
    # shutil.rmtree(os.path.join(out, 'goals'))
    # os.remove(os.path.join(out, 'primers.fasta'))


def test(outdir):
    cli = CLI(".test_out")
    cli.run()
