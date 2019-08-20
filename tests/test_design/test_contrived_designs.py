import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from uuid import uuid4
from pyblast.utils import make_linear, make_circular
from dasi import Design


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

    r1 = random_record(1000)
    r2 = random_record(1000)
    r3 = random_record(1000)

    make_linear([r1, r2, r3])

    goal = r1 + r2 + r3
    make_circular([goal])

    design = Design()
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


# 1. test overlapping fragments
# 2. test overhang pcrs
# 3. test gaps
# 4. 