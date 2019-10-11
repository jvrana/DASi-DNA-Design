from os.path import join

import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import Design


@pytest.mark.parametrize(
    "query", ["pmodkan-ho-pact1-z4-er-vpr.gb", "plko-pef1a-frt-tdtomato-wpre.gb"]
)
def test_num_groups_vs_endpoints(here, paths, query, span_cost):
    primers = make_linear(load_fasta_glob(paths["primers"]))
    templates = load_genbank_glob(paths["templates"])

    query_path = join(here, "data/test_data/genbank/designs", query)
    queries = make_circular(load_genbank_glob(query_path))

    design = Design(span_cost)

    design.add_materials(primers=primers, templates=templates, queries=queries)

    design._blast()
    containers = design.container_list()
    assert len(containers) == 1
    container = containers[0]
    container.expand()
    groups = container.groups()
    print(len(groups) ** 2)

    a_arr = set()
    b_arr = set()

    for g in groups:
        a_arr.add(g.query_region.a)
        b_arr.add(g.query_region.b)

    print(len(a_arr) * len(b_arr))


class TestCompile:

    pass


class TestDesignResult:
    def test_expected_span_length(self, single_processed_results):
        design, results = single_processed_results
        for qk, result in results.items():
            assembly = result.assemblies[0]
            assembly.print()

            assert len(result.query) == sum(assembly.to_df()["span"])

    def test_multi_expected_span_length(self, multi_processed_results):
        design, results = multi_processed_results
        for qk, result in results.items():
            assembly = result.assemblies[0]
            assembly.print()

            assert len(result.query) == sum(assembly.to_df()["span"])


class TestSequenceDesign:
    def test_single_design_sequences(self, single_processed_results):
        design, results = single_processed_results
        for qk, result in results.items():
            result.design_sequences()
            result.design_sequence_output()

    def test_multi_design_sequences(self, multi_processed_results):
        design, results = multi_processed_results
        for qk, result in results.items():
            result.design_sequences()
            result.design_sequence_output()


class TestMultivsSingleProcessedResults:

    pass


class TestMultiProcessing:
    def assembly_equiavlent(self, a1, a2):
        for e1, e2 in zip(a1.edges(), a2.edges()):
            assert e1[0] == e2[0]
            assert e1[1] == e2[1]
            assert e1[2] == e2[2]

    def test_same_results(self, single_processed_results, multi_processed_results):
        design1, results1 = single_processed_results
        design2, results2 = multi_processed_results

        assert len(results1) == len(results2), "Number of results should be the same"

        name_to_key = {design2.seqdb[k].name: k for k in results2}
        key_to_name = {k: design1.seqdb[k].name for k in results1}

        paired_keys = []
        for key1 in results1:
            name1 = key_to_name[key1]
            key2 = name_to_key[name1]
            paired_keys.append((key1, key2))

        for k1, k2 in paired_keys:
            r1 = results1[k1]
            r2 = results2[k2]

            for a1, a2 in zip(r1.assemblies, r2.assemblies):
                self.assembly_equiavlent(a1, a2)
