import dill
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


def test_dump_load_primer_cost(primer_cost):
    dill.loads(dill.dumps(primer_cost))


def test_dump_load_syn_cost(syn_cost):
    dill.loads(dill.dumps(syn_cost))
