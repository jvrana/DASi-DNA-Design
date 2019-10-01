from os.path import join

from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.design import LibraryDesign


def test_library_design(paths, here, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))

    design = LibraryDesign(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile_library()

    # design.optimize_library()

    # design.compile()
    #
    # assert len(design.graphs) == len(queries)
    # assert len(design.graphs) > 1
    #
    # design.optimize()
