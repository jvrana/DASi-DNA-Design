from shoestring.cost import JunctionCost, SynthesisCost
import pytest


@pytest.fixture(scope="module")
def jxncost():
    return JunctionCost()


@pytest.fixture(scope="module")
def syncost(jxncost):
    return SynthesisCost(jxncost)


@pytest.mark.parametrize("span", list(range(-300, 500, 11)))
@pytest.mark.parametrize("ext", [0, 1, 2])
def test_junction_cost(jxncost, span, ext):
    cost = jxncost.junction_cost(span, ext)


def test_plot_junction_cost(jxncost):
    ax = jxncost.plot()


# TODO: do something with design flexibility in cost
def test_plot_flexibility(jxncost):
    jxncost.plot_design_flexibility()


def test_synthesis_cost(syncost):
    syncost.plot()
