import os
import shutil

import pytest

here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture
def outdir(paths):
    out = os.path.join(here, ".out")
    if not os.path.isdir(out):
        os.mkdir(out)
    shutil.copytree(os.path.dirname(paths["templates"]), os.path.join(out, "templates"))
    shutil.copytree(os.path.dirname(paths["registry"]), os.path.join(out, "goals"))
    shutil.copy(paths["primers"], os.path.join(out, "primers.fasta"))
    yield out


def test(outdir):
    print(outdir)
