from shoestring.design import Design
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular


def test_design(paths):
    primers = make_linear(
        load_fasta_glob(paths["primers"])
    )
    templates = load_genbank_glob(paths["templates"])
    queries = make_circular(
        load_genbank_glob(
            paths["queries"]
        )
    )

    design = Design()

    design.add(
        primers=primers,
        templates=templates,
        queries=queries
    )