from shoestring.cost import SpanCost
import pytest
from shoestring.design import Design
from shoestring.cost import SpanCost
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, make_circular
import pytest
from os.path import join


class TestBenchmarkCosts(object):

    def test_span_cost_constructor(self, benchmark):
        benchmark.pedantic(SpanCost, rounds=1, iterations=1)


@pytest.mark.parametrize('query', [
    "pmodkan-ho-pact1-z4-er-vpr.gb",
    'plko-pef1a-frt-tdtomato-wpre.gb'
])
def test_num_groups_vs_endpoints(benchmark, here, paths, query):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, 'data/test_data/genbank/designs', query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design()

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design._blast()

    benchmark(design._blast)