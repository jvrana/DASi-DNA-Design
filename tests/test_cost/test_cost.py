from shoestring.cost import JunctionCost, SynthesisCost


def test_junction_cost():
    JunctionCost()


def test_synthesis_cost():
    s = SynthesisCost(JunctionCost())