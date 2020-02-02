from os.path import join

import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.cost import SpanCost
from dasi.design import Design


@pytest.mark.benchmark
class TestBenchmarkCosts:
    def test_span_cost_constructor(self, benchmark):
        benchmark.pedantic(SpanCost.open, rounds=1, iterations=1)


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize(
    "query", ["pmodkan-ho-pact1-z4-er-vpr.gb", "plko-pef1a-frt-tdtomato-wpre.gb"]
)
def test_benchmark_blast(benchmark, here, paths, query):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design()

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design._blast()

    benchmark(design._blast)
