from shoestring.cost import SpanCost


def test_span_cost_constructor(benchmark):
    benchmark.pedantic(SpanCost, rounds=1, iterations=1)