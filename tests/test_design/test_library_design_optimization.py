import json
from copy import deepcopy
from os.path import join

import networkx as nx
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.constants import Constants
from dasi.design import LibraryDesign
from dasi.models import AlignmentContainer
from dasi.utils import sort_with_keys


def test_library_design_to_df(paths, here, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/library_designs/*.gb")
    queries = make_circular(load_genbank_glob(query_path))
    queries = queries

    design = LibraryDesign(span_cost=span_cost)
    design.n_jobs = 1
    design.add_materials(primers=primers, templates=templates, queries=queries)
    design.compile_library()
    results = design.optimize_library()

    a, b, c = design.to_df()
    a.to_csv("library_design.csv")
    b.to_csv("library_summary.csv")
    with open("designs.json", "w") as f:
        json.dump(c, f)
    print(a)
    print(b)
    print(c)
    # for qk, result in results.items():
    #     df = result.assemblies[0].to_df()
    #     print(design.seqdb[qk].name)
    #     print(df)
