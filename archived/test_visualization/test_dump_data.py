from dasi.design import Design
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
import pytest
from os.path import join
import json


@pytest.mark.parametrize(
    "query",
    [
        # "pmodkan-ho-pact1-z4-er-vpr.gb",
        "plko-pef1a-frt-tdtomato-wpre.gb"
    ],
)
def test_dump_for_visualization(here, paths, query):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design()

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design.compile()

    with open(join(here, "..", "static", "data.json"), "w") as f:
        for container in design.container_list():
            query_region = container.alignment_groups[0].query_region
            data = {
                "query": {
                    "length": query_region.context_length,
                    "name": query_region.name,
                },
                "subjects": [],
            }

            groups = sorted(
                container.groups(),
                key=lambda g: (g.query_region.a, len(g.query_region)),
            )

            for group in groups:
                data["subjects"].append(
                    {
                        "a": group.query_region.a,
                        "b": group.query_region.b,
                        "name": group.alignments[0].subject_region.name,
                    }
                )
            json.dump(data, f)
