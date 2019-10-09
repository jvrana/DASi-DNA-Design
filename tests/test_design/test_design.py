from os.path import join

import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import Design


@pytest.mark.parametrize(
    "query",
    [
        "pins-0a-psv40-citrine-wpre.gb",
        "pmodkan-ho-pact1-z4-er-vpr.gb",
        "plko-pef1a-frt-tdtomato-wpre.gb",
        "pins-01-hu6-sv40-nt1-optgrna.gb",
        "pins-01-hu6-r1-optgrna.gb",
        "pins-01-hu6-r5-optgrna.gb",
    ],
)
def test_real_design(here, paths, query, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"]) + load_genbank_glob(
        paths["registry"]
    )

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost=span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    assert len(design.graphs) == len(queries)
    assert len(design.graphs) == 1

    results = design.optimize()

    for qk, result in results.items():
        result.design_sequences()


# @pytest.mark.parametrize("query", ["plko-pef1a-frt-tdtomato-wpre.gb"])
# def test_real_design(here, paths, query, span_cost):
#     primers = make_linear(load_fasta_glob(paths["primers"]))
#     templates = load_genbank_glob(paths["templates"]) + load_genbank_glob(
#         paths["registry"]
#     )
#
#     query_path = join(here, "data/test_data/genbank/designs", query)
#     queries = make_circular(load_genbank_glob(query_path))
#
#     design = Design(span_cost=span_cost)
#
#     design.add_materials(primers=primers, templates=templates, queries=queries)
#
#     design.compile()
#
#     assert len(design.graphs) == len(queries)
#     assert len(design.graphs) == 1
#
#     results = design.optimize()
#
#     for qk, result in results.items():
#         for a in result.assemblies:
#             print(a.to_df())
#         result.design_sequences()
