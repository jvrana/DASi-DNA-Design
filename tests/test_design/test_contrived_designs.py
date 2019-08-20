import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from uuid import uuid4
from pyblast.utils import make_linear, make_circular
from dasi import Design
from dasi.cost import SpanCost


spancost = SpanCost()


def random_seq(len):
    bases = "AGTC"


    seq = ''
    for _ in range(len):
        i = random.randint(0, 3)
        seq += bases[i]
    return seq


def random_record(len):
    return SeqRecord(Seq(random_seq(len)), id=str(uuid4()))


def test_design_task_1():
    """Fragments with no gaps"""

    r1 = random_record(1000)
    r2 = random_record(1000)
    r3 = random_record(1000)

    make_linear([r1, r2, r3])

    goal = r1 + r2 + r3
    make_circular([goal])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2, r3],
        queries=[goal]
    )
    design.compile()

    for e in list(design.graphs.values())[0].edges(data=True):
        print(e)
    df1, df2 = design.design()
    print(df1)
    assert len(df1) == 6


def test_design_task_2():
    """Fragments with gaps"""
    r1 = random_record(1000)
    r2 = random_record(1000)
    r3 = random_record(1000)

    make_linear([r1, r2, r3])

    goal = r1 + random_record(20)+ r2 + random_record(20) + r3
    make_circular([goal])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2, r3],
        queries=[goal]
    )
    design.compile()

    for e in list(design.graphs.values())[0].edges(data=True):
        print(e)
    df1, df2 = design.design()
    print(df1)
    assert len(df1) == 6


def test_design_task_3():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-100:] + goal[:1000]
    r2 = goal[900:2000]
    r3 = goal[1900:]
    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2, r3],
        queries=[goal]
    )
    design.compile()

    span_lens = []

    for n1, n2, edata in list(design.graphs.values())[0].edges(data=True):
        span_lens.append(edata['span_length'])

    for e in list(design.graphs.values())[0].edges(data=True):
        print(e)
    df1, df2 = design.design()
    d = df1.to_dict()
    print(df1)
    assert len(df1) == 6


def test_design_task4():
    """Test with primers"""
    raise NotImplementedError()