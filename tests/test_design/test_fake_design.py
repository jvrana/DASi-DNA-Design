"""test_fake_design.py.

This performs randomized tests on DASi by generating randomized libraries.

These tests ought to cover a majority of design cases.
"""
import json
from typing import Callable
from typing import Union

import pytest

from dasi import Design
from dasi import LibraryDesign


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


def test_fake_design_short(cached_span_cost):
    design = Design.fake(span_cost=cached_span_cost, n_designs=1)
    run(design)


@parametrize_designs
@pytest.mark.parametrize("design_cls", [LibraryDesign, Design])
def test_fake_design(cached_span_cost, n_designs, synth_prob, synth_size, design_cls):
    design = design_cls.fake(
        span_cost=cached_span_cost,
        n_designs=n_designs,
        n_cyclic_seqs=50,
        n_linear_seqs=50,
        n_primers_from_templates=25,
        random_chunk_prob_int=(synth_prob, synth_prob),
        random_chunk_size_int=(synth_size, synth_size),
    )
    run(design)
