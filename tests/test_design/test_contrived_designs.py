import pytest
import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from uuid import uuid4
from pyblast.utils import make_linear, make_circular
from dasi import Design
from dasi.cost import SpanCost
from itertools import zip_longest
from more_itertools import pairwise

spancost = SpanCost()

# TODO: contrived tests should perform as good or better than expected path

def random_seq(len):
    bases = "AGTC"


    seq = ''
    for _ in range(len):
        i = random.randint(0, 3)
        seq += bases[i]
    return seq


def random_record(len):
    return SeqRecord(Seq(random_seq(len)), id=str(uuid4()))

def print_edge_cost(path, graph):
    total = 0
    path = path[:] + path[:1]
    for n1, n2 in pairwise(path):
        try:
            edata = graph[n1][n2]
            total += edata['weight']
            print((n1, n2, edata['weight']))
        except:
            print((n1, n2, "MISSING EDGE"))
            total += 10**6

    print("TOTAL: {}".format(total))
    return total

def check_design_result(design, expected_path, check_cost=True, check_path=True, path_func=None):
    design.compile()
    path_dict = design.optimize()
    paths = list(path_dict.values())[0]
    best_path = paths[0]

    for x, y in zip_longest(best_path, expected_path):
        z = x == y
        print("{} {} {}".format(z, x, y))

    print("=== BEST PATH COST ===")
    graph = list(design.graphs.values())[0]
    cost1 = print_edge_cost(best_path, graph)

    print("=== EXPECTED PATH COST ===")
    cost2 = print_edge_cost(expected_path, graph)

    print("Num groups: {}".format(len(design.container_list()[0].groups())))

    if check_cost:
        assert cost1 <= cost2

    if check_path:
        if path_func:
            p1 = [path_func(x) for x in best_path]
            p2 = [path_func(x) for x in expected_path]
        else:
            p1 = best_path
            p2 = expected_path
        assert p1 == p2
    return design


def test_blast_has_same_results():
    goal = random_record(3000)
    make_circular([goal])

    r1 = random_record(100) + goal[1000:2000] + random_record(100)
    p1 = goal[970:1030]
    p2 = goal[1970:2030].reverse_complement()
    r2 = goal[2000:] + goal[:1000]
    p3 = goal[1970:2030]

    make_linear([r1, p1, p2, r2, p3])

    size_of_groups = []
    for i in range(20):

        design = Design(spancost)
        design.logger.set_level("INFO")
        design.add_materials(
            primers=[p1, p2, p3],
            templates=[r1, r2],
            queries=[goal],
            fragments=[]
        )

        design.compile()

        for container in design.container_list():
            size_of_groups.append(len(container.groups()))
    assert len(size_of_groups) == 20
    assert len(set(size_of_groups)) == 1


def test_design_with_no_gaps():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[:1000]
    r2 = goal[1000:2000]
    r3 = goal[2000:]
    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2, r3],
        queries=[goal],
        fragments=[],
    )

    expected_path = [
        (0, True, 'A', False),
        (1000, True, 'B', False),
        (1000, True, 'A', False),
        (2000, True, 'B', False),
        (2000, True, 'A', False),
        (3000, True, 'B', False),
    ]

    check_design_result(design, expected_path, path_func=lambda x: (x[0], x[2]))


def test_design_with_overlaps():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1000]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        fragments=[r1, r2, r3],
        queries=[goal],
        templates=[]
    )

    expected_path = [
        (950, False, 'A', True),
        (2000, False, 'B', True),
        (1950, False, 'A', True),
        (3000, False, 'B', True),
        (3000 - 40, False, 'A', True),
        (1000, False, 'B', True),
    ]

    check_design_result(design, expected_path)


def test_design_with_overlaps():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1000]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        fragments=[r1, r2, r3],
        queries=[goal],
        templates=[]
    )

    expected_path = [
        (950, False, 'A', True),
        (2000, False, 'B', True),
        (1950, False, 'A', True),
        (3000, False, 'B', True),
        (3000 - 40, False, 'A', True),
        (1000, False, 'B', True),
    ]

    check_design_result(design, expected_path)


def test_design_with_overlaps_with_templates():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1000]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        fragments=[],
        queries=[goal],
        templates=[r1, r2, r3]
    )

    expected_path = [
        (950, True, 'A', True),
        (2000, True, 'B', True),
        (1950, True, 'A', True),
        (3000, True, 'B', True),
        (3000 - 40, True, 'A', True),
        (1000, True, 'B', True),
    ]

    check_design_result(design, expected_path, check_path=False)


def test_design_task_with_gaps():
    """Fragments with overlaps"""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[:950]
    r2 = goal[1000:2000]
    r3 = goal[2050:]

    make_linear([r1, r2, r3])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2, r3],
        queries=[goal],
        fragments=[]
    )

    expected_path = [
        (0, True, 'A', False),
        (950, True, 'B', False),
        (1000, True, 'A', False),
        (2000, True, 'B', False),
        (2050, True, 'A', False),
        (3000, True, 'B', False),
    ]

    check_design_result(design, expected_path)


# TODO: why are BLAST results so stochastic?
@pytest.mark.parametrize('repeat', range(10))
def test_design_with_overhang_primers(repeat):

    goal = random_record(3000)
    make_circular([goal])

    r1 = random_record(100) + goal[1000:2000] + random_record(100)
    p1 = goal[970:1030]
    p2 = goal[1970:2030].reverse_complement()
    r2 = goal[2000:] + goal[:1000]
    p3 = goal[1970:2030]

    make_linear([r1, p1, p2, r2, p3])

    design = Design(spancost)
    design.add_materials(
        primers=[p1, p2, p3],
        templates=[r1, r2],
        queries=[goal],
        fragments=[]
    )

    expected_path = [
        (970, False, 'A', True),
        (2030, False, 'B', True),
        (1970, False, 'A', True),
        (1000, True, 'B', True),
    ]

    check_design_result(design, expected_path)



def test_requires_synthesis():
    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[1000:2000]
    r2 = goal[200:500]

    make_linear([r1, r2])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2],
        queries=[goal],
        fragments=[]
    )

    expected_path = [
        (200, True, 'A', False),
        (500, True, 'B', False),
        (1000, True, 'A', False),
        (2000, True, 'B', False),
    ]

    check_design_result(design, expected_path)


def test_requires_synthesis_with_template_over_origin():
    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[1000:2000]
    r2 = goal[2500:] + goal[:500]

    make_linear([r1, r2])

    design = Design(spancost)
    design.add_materials(
        primers=[],
        templates=[r1, r2],
        queries=[goal],
        fragments=[]
    )

    expected_path = [
        (500, True, 'B', False),
        (1000, True, 'A', False),
        (2000, True, 'B', False),
        (2500, True, 'A', False),
    ]

    check_design_result(design, expected_path)