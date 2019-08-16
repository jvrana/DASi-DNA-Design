from shoestring.design import LibraryDesign
from shoestring.cost import SpanCost
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
from os.path import join

span_cost = SpanCost()


def test_library_design(paths, here):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, 'data/test_data/genbank/designs/*.gb')
    queries = make_circular(load_genbank_glob(query_path))

    design = LibraryDesign(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    # design.optimize_library()

    # design.compile()
    #
    # assert len(design.graphs) == len(queries)
    # assert len(design.graphs) > 1
    #
    # design.optimize()