"""test_fake_design.py.

This performs randomized tests on DASi by generating randomized libraries.

These tests ought to cover a majority of use cases.
"""
import json
from typing import Callable
from typing import Union

import pytest

from dasi import Design
from dasi import LibraryDesign
from dasi.cost import cached_span_cost as cached_span_cost_default
from dasi.exceptions import DasiDesignException


def parametrize_designs(f: Callable):
    f = pytest.mark.parametrize("n_designs", [1, 3], ids=["1 design", "3 designs"])(f)
    f = pytest.mark.parametrize(
        "synth_prob",
        [0, 0.3, 0.5, 0.9],
        ids=["0% synthetic", "30% synthetic", "50% synthetic", "90% synthetic"],
    )(f)
    f = pytest.mark.parametrize(
        "synth_size", [100, 500, 1000], ids=["100bp", "500bp", "1000bp"]
    )(f)
    return f


def run(design: Union[LibraryDesign, Design]):
    design.compile()
    for qk, graph in design.graphs.items():
        print(graph.number_of_edges())
    design.optimize()
    status = design.status

    print("#" * 25)
    print("# Status")
    print("#" * 25)
    print(json.dumps(status, indent=2))
    print()
    for qk, result in design.results.items():
        for a in result.assemblies:
            print(a)

    react_df, summ_df, design_json = design.to_df()

    print("#" * 25)
    print("### Summary")
    print("#" * 25)
    print(summ_df)
    print()

    print("#" * 25)
    print("### Summary")
    print("#" * 25)
    print(react_df)
    print()


@pytest.mark.parametrize("span_cost", [cached_span_cost_default(), None])
@pytest.mark.parametrize("design_cls", [Design, LibraryDesign])
def test_span_cost_init(span_cost, design_cls):
    """Test if we can use span_cost=None to initialize a new design."""
    design = design_cls.fake(span_cost=span_cost, n_designs=3)
    design.run()


@pytest.mark.parametrize("design_cls", [Design, LibraryDesign])
def test_report(design_cls):
    """Test if we can use span_cost=None to initialize a new design."""
    design = design_cls.fake(n_designs=3)
    design.run()
    design.report().plot_coverage()


@pytest.mark.parametrize("design_cls", [Design, LibraryDesign])
@pytest.mark.parametrize("n_jobs", [1, 3])
@pytest.mark.parametrize("n_designs", [1, 3])
def test_run_with_n_jobs(design_cls, n_jobs, n_designs):
    design = design_cls.fake(n_designs=n_designs)
    design.run(n_jobs=n_jobs)


# TODO: fix library generation for faked plasmids with simliar sequences
def test_cost_comparison_library():
    design1, library = Design.fake(
        n_designs=3,
        n_linear_seqs=50,
        n_cyclic_seqs=50,
        n_primers_from_templates=500,
        shared_length=500,
        return_with_library=True,
    )
    design2 = LibraryDesign(seqdb=design1.seqdb)

    designs = library["design"]
    plasmids = library["cyclic"]
    fragments = library["linear"]
    primers = library["short"]

    design2.add_materials(
        primers=primers, fragments=fragments, templates=plasmids, queries=designs
    )

    design1.run()
    design2.run()

    print("#" * 10 + "\nDesign\n" + "#" * 10)
    print(json.dumps(design1.status, indent=2))

    print("#" * 10 + "\nLibraryDesign\n" + "#" * 10)
    print(json.dumps(design2.status, indent=2))

    print("%" * 10 + "\nDesign Cost\n" + "%" * 10)
    for qk, s in design1.status.items():
        print(s["assemblies"])

    print("%" * 10 + "\nLibraryDesign Cost\n" + "%" * 10)
    for qk, s in design2.status.items():
        print(s["assemblies"])

    design2.report().plot_coverage(show=True)
    print(design2.to_df()[1])


@pytest.mark.parametrize("n_designs", [1, 3])
@pytest.mark.parametrize("design_cls", [Design, LibraryDesign])
def test_design_status(cached_span_cost, n_designs, design_cls):
    """Tests design status reporting."""
    design = design_cls.fake(span_cost=cached_span_cost, n_designs=n_designs)

    assert not design.query_keys
    design.compile()
    assert design.query_keys
    for qk in design.query_keys:
        assert design.status[qk]["run"] is False
        assert design.status[qk]["compiled"] is True

    design.optimize()
    for qk in design.query_keys:
        assert design.status[qk]["run"] is True
        assert design.status[qk]["compiled"] is True


def test_design_optimize_cannot_run_before_compile(cached_span_cost):
    design = Design.fake(n_designs=1)
    with pytest.raises(DasiDesignException):
        design.optimize()
    design.compile()
    design.optimize()


@pytest.mark.parametrize(
    "synth_prob",
    [0, 0.3, 0.5, 0.9],
    ids=["0% synthetic", "30% synthetic", "50% synthetic", "90% synthetic"],
)
@pytest.mark.parametrize("n_designs", [1, 3], ids=["1 design", "3 designs"])
@pytest.mark.parametrize(
    "synth_size", [100, 500, 1000], ids=["100bp", "500bp", "1000bp"]
)
@pytest.mark.parametrize("shared_length", [0, 1000])
@pytest.mark.parametrize("design_cls", [LibraryDesign, Design])
def test_fake_design(
    cached_span_cost, n_designs, synth_prob, synth_size, design_cls, shared_length
):
    design = design_cls.fake(
        span_cost=cached_span_cost,
        n_designs=n_designs,
        n_cyclic_seqs=50,
        n_linear_seqs=50,
        n_primers_from_templates=500,
        random_chunk_prob_int=(synth_prob, synth_prob),
        random_chunk_size_int=(synth_size, synth_size),
        design_sequence_similarity_length=shared_length,
    )
    run(design)
