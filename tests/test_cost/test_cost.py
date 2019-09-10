from dasi.cost import (
    PrimerCostBuilder,
    SynthesisCostBuilder,
    PrimerParams,
    SynthesisParams,
)
import pytest
import pylab as plt
import numpy as np


@pytest.fixture(scope="module")
def primer_cost():
    return PrimerCostBuilder.from_params(PrimerParams)


@pytest.fixture(scope="module")
def syn_cost(primer_cost):
    return SynthesisCostBuilder.from_params(SynthesisParams, primer_cost)


class TestPlotters(object):
    def test_plot_primer_cost(self, primer_cost):
        primer_cost.plot()
        plt.show()

    def test_plot_syn_cost(self, syn_cost):
        syn_cost.plot()
        plt.show()


def test_primer_cost(primer_cost):
    primer_cost(100, (1, 1))


def test_primer_cost_df(primer_cost):
    primer_cost.to_df()


def test_syn_cost(syn_cost):
    syn_cost(100, (1, 1))


def test_syn_cost_df(syn_cost):
    df = syn_cost.to_df()


def test_syn_cost_retrieve_many_times(syn_cost):
    syn_cost(np.arange(1000), (1, 1))


# @pytest.fixture(scope="module")
# def jxncost():
#     return JunctionCost()
#
#
# @pytest.fixture(scope="module")
# def syncost(jxncost):
#     return SynthesisCost(jxncost)
#
#
# @pytest.mark.parametrize("span", list(range(-300, 500, 11)))
# @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (1, 1)])
# def test_junction_cost(jxncost, span, ext):
#     cost = jxncost.cost(span, ext)
#
#
# def test_plot_junction_cost(jxncost):
#     ax = jxncost.plot()
#
#
# # TODO: do something with design flexibility in cost
# def test_plot_flexibility(jxncost):
#     jxncost.plot_design_flexibility()
#
#
# def test_synthesis_cost(syncost):
#     syncost.plot()
#
#
# def test_span_cost():
#     span_cost = SpanCost()
#     span_cost.plot()
# #     x = np.arange(-500, 3000)
# #     y = span_cost.cost(x, (1, 1))
# #
# #     sns.lineplot(x=x, y=y)
# #     plt.show()
