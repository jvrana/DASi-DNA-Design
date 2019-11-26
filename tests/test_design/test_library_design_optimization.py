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
    queries = [queries[-1]]

    design = LibraryDesign(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile_library()

    results = design.optimize_library()

    result = list(results.values())[0]
    df = result.assemblies[0].to_df()

    print(df)
