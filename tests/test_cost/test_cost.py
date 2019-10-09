import os

import numpy as np
import pylab as plt
import pytest

from dasi.cost import PrimerCostModel
from dasi.cost import PrimerParams
from dasi.cost import SpanCost
from dasi.cost import SynthesisCostModel
from dasi.cost import SynthesisParams


@pytest.fixture(scope="module")
def primer_cost():
    return PrimerCostModel.from_params(PrimerParams)


@pytest.fixture(scope="module")
def syn_cost(primer_cost):
    return SynthesisCostModel.from_params(SynthesisParams, primer_cost)


@pytest.fixture(scope="module")
def span_cost(primer_cost, syn_cost):
    return SpanCost(syn_cost)


class TestPlotters:
    def test_plot_primer_cost(self, primer_cost):
        primer_cost.plot()
        plt.show()

    def test_plot_syn_cost(self, syn_cost):
        syn_cost.plot()
        plt.show()

    def test_plot_span_cost(self, span_cost):
        span_cost.plot()
        plt.show()


class TestDf:
    def test_primer_cost_df(self, primer_cost):
        primer_cost.to_df()

    def test_syn_cost_df(self, syn_cost):
        syn_cost.to_df()

    def test_span_cost_df(self, span_cost):
        span_cost.to_df()


class TestCall:
    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_primer_cost_df(self, primer_cost, ext):
        primer_cost(np.arange(-300, 1000), ext)

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_syn_cost_df(self, syn_cost, ext):
        syn_cost(np.arange(-300, 1000), (0, 0))

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_span_cost_df(self, span_cost, ext):
        span_cost(np.arange(-300, 1000), (0, 0))


class TestEdgeCases:
    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_primer_extreme_left(self, primer_cost, ext):
        df = primer_cost(-1000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_primer_extreme_right(self, primer_cost, ext):
        df = primer_cost(1000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_syn_extreme_left(self, syn_cost, ext):
        df = syn_cost(-1000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_syn_extreme_right(self, syn_cost, ext):
        df = syn_cost(5000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_span_extreme_left(self, span_cost, ext):
        df = span_cost.cost(-1000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf

    @pytest.mark.parametrize("ext", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_span_extreme_right(self, span_cost, ext):
        df = span_cost.cost(5000, ext)
        print(df.col["cost"])
        assert df.data["cost"][0] == np.inf


class TestSerialization:

    here = os.path.abspath(os.path.dirname(__file__))

    def test_dumpb(self, span_cost):
        span_cost.dumpb()

    def test_loadb(self, span_cost):
        s = span_cost.dumpb()
        span_cost2 = SpanCost.loadb(s)
        assert len(span_cost2.to_df()) == len(span_cost.to_df())

    def test_dump_and_load(self, tmp_path, span_cost):
        tmp_file = os.path.join(tmp_path, "span_cost.b")
        span_cost.dump(tmp_file)
        assert os.path.isfile(tmp_file)
        span_cost2 = SpanCost.load(tmp_file)
        assert len(span_cost2.to_df()) == len(span_cost.to_df())
