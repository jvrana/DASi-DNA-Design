from shoestring.cost import JunctionCost, SynthesisCost, SpanCost
import pytest
import numpy as np
import seaborn as sns
import pylab as plt


@pytest.fixture(scope="module")
def jxncost():
    return JunctionCost()


@pytest.fixture(scope="module")
def syncost(jxncost):
    return SynthesisCost(jxncost)


@pytest.mark.parametrize("span", list(range(-300, 500, 11)))
@pytest.mark.parametrize("ext", [(0, 0), (1, 0), (1, 1)])
def test_junction_cost(jxncost, span, ext):
    cost = jxncost.cost(span, ext)


def test_plot_junction_cost(jxncost):
    ax = jxncost.plot()


# TODO: do something with design flexibility in cost
def test_plot_flexibility(jxncost):
    jxncost.plot_design_flexibility()


def test_synthesis_cost(syncost):
    syncost.plot()


def test_span_cost():
    span_cost = SpanCost()
    x = np.arange(-500, 3000)
    y = span_cost.cost(x, (1, 1))

    sns.lineplot(x=x, y=y)
    plt.show()
