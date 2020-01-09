import shutil
from os import mkdir
from os import remove
from os.path import abspath
from os.path import dirname
from os.path import isdir
from os.path import isfile
from os.path import join

import pytest

from dasi.command_line import DasiCLI

here = abspath(dirname(__file__))


@pytest.fixture
def outdir(paths):
    out = join(here, ".test_out")
    if not isdir(out):
        mkdir(out)
    templates = join(out, "templates")
    goals = join(out, "goals")
    primers = join(out, "primers.fasta")

    if isdir(templates):
        shutil.rmtree(templates)

    if isdir(goals):
        shutil.rmtree(goals)

    if isfile(primers):
        remove(primers)

    shutil.copytree(dirname(paths["registry"]), templates)
    shutil.copytree(dirname(paths["queries"]), goals)
    shutil.copy(paths["primers"], join(out, primers))
    yield out
    # shutil.rmtree(join(out, 'templates'))
    # shutil.rmtree(join(out, 'goals'))
    # os.remove(join(out, 'primers.fasta'))


def test_cli(outdir):
    cli = DasiCLI(join(here, ".test_out"))
    cli.run()
