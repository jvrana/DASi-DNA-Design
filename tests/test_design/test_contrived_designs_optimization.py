import itertools
import random
from typing import List
from typing import Tuple
from uuid import uuid4

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import Design
from dasi.cost import SpanCost
from dasi.design.graph_builder import AssemblyNode


def random_seq(n_bases):
    bases = "AGTC"
    seq = ""
    for _ in range(n_bases):
        i = random.randint(0, 3)
        seq += bases[i]
    return seq


def random_record(n_bases):
    return SeqRecord(Seq(random_seq(n_bases)), id=str(uuid4()))


def print_edge_cost(path, graph):
    total = 0
    path = path[:] + path[:1]
    for n1, n2 in pairwise(path):
        try:
            edata = graph[n1][n2]
            total += edata["weight"]
            print((n1, n2, edata))
        except KeyError:
            print((n1, n2, "MISSING EDGE"))
            total += np.inf

    print("TOTAL: {}".format(total))
    return total


def change_index(n, i):
    return (i, n[1], n[2], n[3])


def fix_path(path, graph, length):
    for n1, n2 in pairwise(path):
        edata = graph.get_edge_data(n1, n2)
        if edata is None:
            n1_arr = []
            n2_arr = []
            edata_arr = []
            if n1.index <= length:
                n1_arr.append(change_index(n1, n1.index - length))
            if n1.index >= length:
                n1_arr.append(change_index(n1, n1.index - length))

            if n2.index <= length:
                n2_arr.append(change_index(n2, n2.index - length))
            if n2.index >= length:
                n2_arr.append(change_index(n2, n2.index - length))

            for _n1, _n2 in itertools.product(n1_arr, n2_arr):
                _edata = graph.get_edge_data(_n1, _n2)
                if _edata:
                    edata_arr.append((_n1, _n2, _edata))

            edata_arr = sorted(edata_arr, key=lambda x: x[2]["weight"])

            if edata_arr:
                return edata_arr[0][0], edata_arr[0][1]


class NoSolution(Exception):
    pass


def check_design_result(
    design, expected_path: List[Tuple], check_cost=False, check_path=True
):

    expected_path = [AssemblyNode(*t) for t in expected_path]

    # compile the design
    design.compile()

    # get the results
    results = list(design.optimize().values())
    result = results[0]
    assemblies = result.assemblies

    if not assemblies:
        raise NoSolution("There are no solution.")

    best_solution = assemblies[0]
    expected_solution = result._add_assembly_from_path(expected_path)

    dist, explain = best_solution.edit_distance(expected_solution, explain=True)

    print("Best nodes")
    for n in best_solution._nodes:
        print(n)

    print("=== DIFF ===")
    best_solution.print_diff(expected_solution)

    print("=== BEST ===")
    df1 = best_solution.to_df()
    print(df1)

    print("=== EXPECTED ===")
    df2 = expected_solution.to_df()
    print(df2)

    print("Best: {}".format(best_solution.cost()))
    print("Expected: {}".format(expected_solution.cost()))

    if check_path:
        # check DataFrame
        cols = ["query_start", "query_end"]
        assert df1[[*cols]].equals(df2[[*cols]])

        # check edges
        for e1, e2 in zip(
            best_solution.edges(data=False), expected_solution.edges(data=False)
        ):
            assert e1 == e2, "{} != {}".format(e1, e2)

    if check_cost:
        assert best_solution.cost() <= expected_solution.cost()
        assert expected_solution.cost() != np.inf


def test_blast_has_same_results(span_cost):
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

        design = Design(span_cost)
        design.logger.set_level("INFO")
        design.add_materials(
            primers=[p1, p2, p3], templates=[r1, r2], queries=[goal], fragments=[]
        )

        design.compile()

        for container in design.container_list():
            size_of_groups.append(len(container.groups()))
    assert len(size_of_groups) == 20
    assert len(set(size_of_groups)) == 1


def test_design_with_no_gaps(span_cost):
    """Fragments with overlaps."""

    goal = random_record(3000)
    make_circular([goal])

    r1 = SeqRecord(Seq("NNNNNN")) + goal[:1000] + SeqRecord(Seq("NNNNNN"))
    r2 = goal[1000:2000]
    r3 = goal[2000:]
    make_linear([r1, r2, r3])

    design = Design(span_cost)
    design.add_materials(
        primers=[], templates=[r1, r2, r3], queries=[goal], fragments=[]
    )

    expected_path = [
        (0, True, "A", False),
        (1000, True, "B", False),
        (1000, True, "A", False),
        (2000, True, "B", False),
        (2000, True, "A", False),
        (3000, True, "B", False),
    ]

    check_design_result(design, expected_path)


def test_design_with_overlaps(span_cost):
    """Fragments with overlaps."""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1000]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(span_cost)
    design.add_materials(
        primers=[], fragments=[r1, r2, r3], queries=[goal], templates=[]
    )

    expected_path = [
        (950, False, "A", True),
        (2000, False, "B", True),
        (1950, False, "A", True),
        (3000, False, "B", True),
        (3000 - 40, False, "A", True),
        (1000, False, "B", True),
    ]

    check_design_result(design, expected_path)


def test_design_with_overlaps2(span_cost):
    """Fragments with overlaps."""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1001]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(span_cost)
    design.add_materials(
        primers=[], fragments=[r1, r2, r3], queries=[goal], templates=[]
    )

    expected_path = [
        (950, False, "A", True),
        (2000, False, "B", True),
        (1950, False, "A", True),
        (3000, False, "B", True),
        (3000 - 40, False, "A", True),
        (4001, False, "B", True),
    ]

    check_design_result(design, expected_path)


def test_design_with_overlaps_with_templates(span_cost):
    """Fragments with overlaps."""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[-40:] + goal[:1000]
    r2 = goal[950:2000]
    r3 = goal[1950:]
    make_linear([r1, r2, r3])

    design = Design(span_cost)
    design.add_materials(
        primers=[], fragments=[], queries=[goal], templates=[r1, r2, r3]
    )

    expected_path = [
        (950, True, "A", True),
        (2000, True, "B", True),
        (1950, True, "A", True),
        (3000, True, "B", True),
        (3000 - 40, True, "A", True),
        (4000, True, "B", True),
    ]

    check_design_result(design, expected_path, check_path=False)


def test_design_task_with_gaps(span_cost):
    """Fragments with overlaps."""

    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[:950]
    r2 = goal[1000:2000]
    r3 = goal[2050:]

    make_linear([r1, r2, r3])

    design = Design(span_cost)
    design.add_materials(
        primers=[], templates=[r1, r2, r3], queries=[goal], fragments=[]
    )

    expected_path = [
        (0, True, "A", False),
        (950, True, "B", False),
        (1000, True, "A", False),
        (2000, True, "B", False),
        (2050, True, "A", False),
        (3000, True, "B", False),
    ]

    check_design_result(design, expected_path)


@pytest.mark.parametrize("repeat", range(3))
def test_design_with_overhang_primers(repeat, span_cost):

    goal = random_record(3000)
    make_circular([goal])

    r1 = random_record(100) + goal[1000:2000] + random_record(100)
    p1 = goal[970:1030]
    p2 = goal[1970:2030].reverse_complement()
    r2 = goal[2000:] + goal[:1000]
    p3 = goal[1970:2030]

    make_linear([r1, p1, p2, r2, p3])

    design = Design(span_cost)
    design.add_materials(
        primers=[p1, p2, p3], templates=[r1, r2], queries=[goal], fragments=[]
    )

    expected_path = [
        (970, False, "A", True),
        (2030, False, "B", True),
        (1970, False, "A", True),
        (4000, True, "B", True),
    ]

    check_design_result(design, expected_path)


def test_amplify_fragment(span_cost):
    goal = random_record(4000)
    make_circular([goal])

    f1 = goal[500:3000]
    f2 = goal[1000:4000]
    f3 = goal[:500]


def test_requires_synthesis(span_cost):
    goal = random_record(4000)
    make_circular([goal])

    r1 = goal[1000:2000]
    r2 = goal[200:500]

    make_linear([r1, r2])

    design = Design(span_cost)
    design.add_materials(primers=[], templates=[r1, r2], queries=[goal], fragments=[])

    expected_path = [
        (200, True, "A", False),
        (500, True, "B", False),
        (1000, True, "A", False),
        (2000, True, "B", False),
    ]

    check_design_result(design, expected_path)


def test_requires_synthesis_with_template_over_origin(span_cost):
    goal = random_record(5000)
    make_circular([goal])

    r1 = goal[1000:2000]
    r2 = goal[3500:] + goal[:500]

    make_linear([r1, r2])

    design = Design(span_cost)
    design.add_materials(primers=[], templates=[r1, r2], queries=[goal], fragments=[])

    expected_path = [
        (500, True, "B", False),
        (1000, True, "A", False),
        (2000, True, "B", False),
        (3500, True, "A", False),
    ]

    check_design_result(design, expected_path)


def test_very_long_synthesizable_region(span_cost):
    goal = random_record(10000)
    make_circular([goal])

    r1 = goal[4177:4255]
    r2 = goal[4188:4225]

    make_linear([r1])

    design = Design(span_cost)
    design.add_materials(primers=[], templates=[r1], queries=[goal], fragments=[])

    expected_path = [
        (500, True, "B", False),
        (1000, True, "A", False),
        (2000, True, "B", False),
        (2500, True, "A", False),
    ]

    with pytest.raises(NoSolution):
        check_design_result(design, expected_path)


def test_single_fragment(span_cost):
    goal = random_record(3000)
    make_circular([goal])

    r1 = goal[177:2255]

    make_linear([r1])

    design = Design(span_cost)
    design.add_materials(primers=[], templates=[r1], queries=[goal], fragments=[])

    expected_path = [(177, True, "A", False), (2255, True, "B", False)]

    check_design_result(design, expected_path)


def test_fully_overlapped(span_cost):
    goal = random_record(2000)
    make_circular([goal])

    r1 = goal[1100:1300]
    p1 = goal[1177 : 1177 + 30]
    p2 = goal[1188 : 1188 + 30]
    p3 = goal[1225 - 30 : 1225].reverse_complement()

    make_linear([r1, p1, p2, p3])

    design = Design(span_cost)
    design.add_materials(
        primers=[p1, p2, p3], templates=[r1], queries=[goal], fragments=[]
    )

    expected_path = [(1177, False, "A", False), (1300, True, "B", False)]

    check_design_result(design, expected_path)


def test_case(span_cost):
    """This is a test case which has previously failed to find a solution.

    The case is that there are just two small fragments with a small
    <10bp gap. The solution should be to PCR amplify both fragment and
    synthesize the rest of the plasmid.
    """
    goal = random_record(2000)
    make_circular([goal])

    r1 = goal[1188:1230]
    r2 = goal[1238:1282]

    make_linear([r1, r2])

    design = Design(span_cost)
    design.add_materials(primers=[], templates=[r1, r2], queries=[goal], fragments=[])

    expected_path = [(1238, True, "A", False), (1282, True, "B", False)]

    check_design_result(design, expected_path)


def test_a_reverse_pcr_fragment(span_cost):

    goal = random_record(3000)
    make_circular([goal])

    t1 = goal[1000:2500].reverse_complement()
    p1 = goal[2500 - 20 : 2510].reverse_complement()

    make_linear([p1, t1])

    design = Design(span_cost)
    design.add_materials(primers=[p1], templates=[t1], queries=[goal], fragments=[])

    expected_path = [(1000, True, "A", False), (2510, False, "B", False)]

    check_design_result(design, expected_path)


def test_case2(span_cost):

    goal = random_record(5000)
    make_circular([goal])

    p1 = goal[1096 : 1096 + 20]
    p2 = goal[694 - 20 : 694].reverse_complement()

    p3 = goal[665 - 20 : 665].reverse_complement()

    t1 = random_seq(100) + goal[1000:] + goal[:700] + random_seq(100)

    make_linear([p1, p2, p3])
    make_circular([t1])

    design = Design(span_cost)

    design.add_materials(
        primers=[p1, p2, p3], templates=[t1], queries=[goal], fragments=[]
    )

    expected_path = [(1238, True, "A", False), (1282, True, "B", False)]

    # check_design_result(design, expected_path)

    design.compile()
    results = design.optimize()
    result = list(results.values())[0]
    for a in result.assemblies:
        print(a)
        print(a.cost())
        for n1, n2, edata in a.edges():
            print(n1, n2)
        # print(a.to_df())
