import pytest

from dasi import Design
from dasi import LibraryDesign


@pytest.mark.parametrize("n_designs", [1, 3])
def test_fake_design(n_designs):
    design = Design.fake(n_designs=n_designs, n_primers_from_templates=10)
    design.compile()
    design.optimize()


@pytest.mark.parametrize("n_designs", [1, 3])
def test_fake_library_design(n_designs):
    design = LibraryDesign.fake(n_designs=n_designs, n_primers_from_templates=10)
    design.compile()
    design.optimize()
