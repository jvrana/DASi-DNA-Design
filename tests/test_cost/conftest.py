import pytest

from dasi import cost


@pytest.fixture(scope="module")
def primer_cost():
    return cost.PrimerCostModel.open()


@pytest.fixture(scope="module")
def syn_cost(primer_cost):
    return cost.SynthesisCostModel.open(primer_cost=primer_cost)


@pytest.fixture(scope="module")
def span_cost(primer_cost, syn_cost):
    return cost.SpanCost(syn_cost)
