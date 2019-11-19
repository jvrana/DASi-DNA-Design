"""Tests that check the values of the span cost."""
import numpy as np
import pytest

from dasi import cost


class TestPrimerCostValues:
    def test_no_primer_cost(self, primer_cost):
        df = primer_cost.cost(np.arange(-20, 100), (0, 0))
        assert np.all(df.data["left_ext"] == 0)
        assert np.all(df.data["material"] == 0)
        assert np.all(df.data["right_ext"] == 0)

    def test_lprimer_cost(self, primer_cost):
        df = primer_cost.cost(np.arange(-20, 100), (1, 0))
        print(df)
        assert np.all(df.data["material"] > 0)
        assert np.all(df.data["right_ext"] == 0)

    def test_rprimer_cost(self, primer_cost):
        df = primer_cost.cost(np.arange(-20, 100), (0, 1))
        assert np.all(df.data["material"] > 0)
        assert np.all(df.data["left_ext"] == 0)

    def test_lrprimer_cost(self, primer_cost):
        df = primer_cost.cost(np.arange(-20, 100), (1, 1))
        print(df)
        assert np.all(df.data["material"] > 0)

    @pytest.mark.parametrize(
        "ext", [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1))]
    )
    def test_relative_primer_cost(self, primer_cost, ext):
        df1 = primer_cost.cost(np.arange(-20, 20), ext[0])
        df2 = primer_cost.cost(np.arange(-20, 20), ext[1])
        print(df1.data["material"])
        print(df2.data["material"])
        a = df1.data["material"] < df2.data["material"]
        print(a)
        assert np.all(a)

    def test_symmetrical_primer_cost(self, primer_cost):
        df1 = primer_cost.cost(np.arange(-20, 100), (1, 0))
        df2 = primer_cost.cost(np.arange(-20, 100), (0, 1))
        assert np.all(df1.data["material"] == df2.data["material"])


class TestBasicPrimerMaterialCost:
    @pytest.fixture(scope="module")
    def small_span_cost(self):
        params = cost.open_params()
        span = np.arange(510, 515)
        span_cost = cost.SpanCost.from_json(params, override_span=span)
        return span_cost

    small_span = pytest.mark.parametrize("span", [-30, -20, np.arange(-30, 100)])

    large_span = pytest.mark.parametrize("span", [512])

    @small_span
    def test_basic_lprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.primer_cost.cost(span, (1, 0))
        print(df)
        assert np.all(df.data["material"] > 0)

    @small_span
    def test_basic_rprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.primer_cost.cost(span, (0, 1))
        print(df)
        assert np.all(df.data["material"] > 0)

    @small_span
    def test_basic_lrprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.primer_cost.cost(span, (1, 1))
        print(df)
        assert np.all(df.data["material"] > 0)

    @small_span
    def test_basic_no_primer_material_cost(self, small_span_cost, span):
        df = small_span_cost.primer_cost.cost(span, (0, 0))
        print(df)
        assert np.all(df.data["material"] == 0)

    @large_span
    def test_syn_lprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.syn_cost.cost(span, (1, 0))
        assert np.all(df.data["lprimer_material"] > 0)
        assert np.all(df.data["rprimer_material"] == 0)

    @large_span
    def test_syn_rprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.syn_cost.cost(span, (0, 1))
        assert np.all(df.data["lprimer_material"] == 0)
        assert np.all(df.data["rprimer_material"] > 0)

    @large_span
    def test_syn_lrprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.syn_cost.cost(span, (1, 1))
        assert np.all(df.data["lprimer_material"] > 0)
        assert np.all(df.data["rprimer_material"] > 0)

    @large_span
    def test_syn_no_primer_material_cost(self, small_span_cost, span):
        df = small_span_cost.syn_cost.cost(span, (0, 0))
        assert np.all(df.data["lprimer_material"] == 0)
        assert np.all(df.data["rprimer_material"] == 0)

    @large_span
    def test_span_lprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.cost(span, (1, 0))
        assert np.all(df.data["lprimer_material"] > 0)
        assert np.all(df.data["rprimer_material"] == 0)

    @large_span
    def test_span_rprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.cost(span, (0, 1))
        assert np.all(df.data["lprimer_material"] == 0)
        assert np.all(df.data["rprimer_material"] > 0)

    @large_span
    def test_span_lrprimer_material_cost(self, small_span_cost, span):
        df = small_span_cost.cost(span, (1, 1))
        assert np.all(df.data["lprimer_material"] > 0)
        assert np.all(df.data["rprimer_material"] > 0)

    @large_span
    def test_span_no_primer_material_cost(self, small_span_cost, span):
        df = small_span_cost.cost(span, (0, 0))
        assert np.all(df.data["lprimer_material"] == 0)
        assert np.all(df.data["rprimer_material"] == 0)


class TestCheckValuesForSpanCost:

    small_span = pytest.mark.parametrize(
        "span",
        [
            (np.arange(500, 505), True),
            (np.arange(100, 1000), True),
            (np.arange(100, 200), True),
            (np.arange(900, 2000), True),
            (np.arange(50, 200), True),
            (np.arange(-100, 100), False),
        ],
        ids=[
            "500,505",
            "100,1000",
            "0,200",
            "900,2000",
            "50,200",
            "-100,100,check_ext",
        ],
    )
    overlapping_span = pytest.mark.parametrize("span", [(-100, 100)])

    @small_span
    def test_no_primer_material_cost(self, cached_span_cost, span):
        span, check_ext = span
        df1 = cached_span_cost.cost(span, (0, 0))
        assert np.all(df1.data["lprimer_material"] == 0)
        assert np.all(df1.data["rprimer_material"] == 0)

    @small_span
    def test_lprimer_material_cost(self, cached_span_cost, span):
        span, check_ext = span
        df1 = cached_span_cost.cost(span, (1, 0))

        a = df1.data["lprimer_material"] > 0
        b = df1.data["rprimer_material"] == 0
        c = df1.data["gene_cost"] > 0
        d = df1.data["left_ext"] > 0
        e = df1.data["right_ext"] == 0

        assert np.all(np.logical_or(~c, a)), (
            "Expect lprimer material to be non-zero " "for values where gene_cost > 0"
        )

        assert np.all(np.logical_or(~c, b)), (
            "Expect rprimer material to be zero " "for values where gene_cost > 0"
        )

        if check_ext:
            assert np.all(np.logical_or(c, d)), (
                "Expect left_ext to be non-zero " "for values where gene_cost == 0"
            )

            assert np.all(np.logical_or(c, e)), (
                "Expect right_ext to be zero " "for values where gene_cost == 0"
            )

    @small_span
    def test_rprimer_material_cost(self, cached_span_cost, span):
        span, check_ext = span
        df1 = cached_span_cost.cost(span, (0, 1))

        a = df1.data["lprimer_material"] == 0
        b = df1.data["rprimer_material"] > 0
        c = df1.data["gene_cost"] > 0
        d = df1.data["left_ext"] == 0
        e = df1.data["right_ext"] > 0

        assert np.all(np.logical_or(~c, a)), (
            "Expect lprimer material to be non-zero " "for values where gene_cost > 0"
        )

        assert np.all(np.logical_or(~c, b)), (
            "Expect rprimer material to be zero " "for values where gene_cost > 0"
        )

        if check_ext:
            assert np.all(np.logical_or(c, d)), (
                "Expect left_ext to be non-zero " "for values where gene_cost == 0"
            )

            assert np.all(np.logical_or(c, e)), (
                "Expect right_ext to be zero " "for values where gene_cost == 0"
            )

    @small_span
    def test_lrprimer_material_cost(self, cached_span_cost, span):
        span, check_ext = span
        df1 = cached_span_cost.cost(span, (1, 1))

        a = df1.data["lprimer_material"] > 0
        b = df1.data["rprimer_material"] > 0
        c = df1.data["gene_cost"] > 0
        d = df1.data["left_ext"] > 0
        e = df1.data["right_ext"] > 0

        assert np.all(np.logical_or(~c, a)), (
            "Expect lprimer material to be non-zero " "for values where gene_cost > 0"
        )

        assert np.all(np.logical_or(~c, b)), (
            "Expect rprimer material to be zero " "for values where gene_cost > 0"
        )

        if check_ext:
            assert np.all(np.logical_or(c, d)), (
                "Expect left_ext to be non-zero " "for values where gene_cost == 0"
            )

            assert np.all(np.logical_or(c, e)), (
                "Expect right_ext to be zero " "for values where gene_cost == 0"
            )

    @small_span
    def test_no_primer_material_cost(self, cached_span_cost, span):
        span, check_ext = span
        df1 = cached_span_cost.cost(span, (0, 0))

        a = df1.data["lprimer_material"] == 0
        b = df1.data["rprimer_material"] == 0
        c = df1.data["gene_cost"] > 0
        d = df1.data["left_ext"] == 0
        e = df1.data["right_ext"] == 0

        assert np.all(np.logical_or(~c, a)), (
            "Expect lprimer material to be non-zero " "for values where gene_cost > 0"
        )

        assert np.all(np.logical_or(~c, b)), (
            "Expect rprimer material to be zero " "for values where gene_cost > 0"
        )

        if check_ext:
            assert np.all(np.logical_or(c, d)), (
                "Expect left_ext to be non-zero " "for values where gene_cost == 0"
            )

            assert np.all(np.logical_or(c, e)), (
                "Expect right_ext to be zero " "for values where gene_cost == 0"
            )
