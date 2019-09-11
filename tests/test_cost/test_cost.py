from dasi.cost import (
    PrimerCostBuilder,
    SynthesisCostBuilder,
    PrimerParams,
    SynthesisParams,
    SpanCost
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

@pytest.fixture(scope="module")
def span_cost(primer_cost, syn_cost):
    return SpanCost(syn_cost)


class TestPlotters(object):
    def test_plot_primer_cost(self, primer_cost):
        primer_cost.plot()
        plt.show()

    def test_plot_syn_cost(self, syn_cost):
        syn_cost.plot()
        plt.show()

    def test_plot_span_cost(self, syn_cost):
        span_cost = SpanCost(syn_cost, syn_cost.primer_cost)
        span_cost.plot()
        plt.show()


class TestDf(object):

    def test_primer_cost_df(self, primer_cost):
        primer_cost.to_df()

    def test_syn_cost_df(self, syn_cost):
        syn_cost.to_df()

    def test_span_cost_df(self, span_cost):
        span_cost.to_df()


class TestCall(object):

    @pytest.mark.parametrize('ext', [
        (0, 0), (1, 0), (0, 1), (1, 1)
    ])
    def test_primer_cost_df(self, primer_cost, ext):
        primer_cost(np.arange(-300, 1000), ext)

    @pytest.mark.parametrize('ext', [
        (0, 0), (1, 0), (0, 1), (1, 1)
    ])
    def test_syn_cost_df(self, syn_cost, ext):
        syn_cost(np.arange(-300, 1000), (0, 0))

    @pytest.mark.parametrize('ext', [
        (0, 0), (1, 0), (0, 1), (1, 1)
    ])
    def test_span_cost_df(self, span_cost, ext):
        span_cost(np.arange(-300, 1000), (0, 0))
