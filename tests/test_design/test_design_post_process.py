import random
from uuid import uuid4

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi import Design
from dasi.models import AssemblyNode
from dasi.utils import Region


def make_linear_and_id(rlist):
    make_linear(rlist)
    for r in rlist:
        r.id = str(uuid4())


def make_circular_and_id(rlist):
    make_circular(rlist)
    for r in rlist:
        r.id = str(uuid4())


def random_seq(n_bases):
    bases = "AGTC"
    seq = ""
    for _ in range(n_bases):
        i = random.randint(0, 3)
        seq += bases[i]
    return seq


def random_record(n_bases):
    return SeqRecord(Seq(random_seq(n_bases)), id=str(uuid4()))


class TestPostProcess:
    def test_add_special_partition_node(self, span_cost):
        """This test adds a new unique node 'n3' with a unique type to simulate
        adding a partitioning to the graph."""
        goal = random_record(4000)
        make_circular_and_id([goal])

        r1 = goal[1000:2000]
        r2 = goal[200:500]

        make_linear_and_id([r1, r2])

        design = Design(span_cost)
        design.add_materials(
            primers=[], templates=[r1, r2], queries=[goal], fragments=[]
        )

        design.compile()

        import networkx as nx

        for qk, g in design.graphs.items():
            query = design.seqdb[qk]
            gcopy = nx.DiGraph(g)
            for n1, n2, edata in g.edges(data=True):
                r = Region(n1.index, n2.index, len(query.seq), cyclic=True)
                if n1.type == "B" and n2.type == "A":
                    # index = int((n1.index + n2.index) / 2)
                    delta = int(len(r) / 2)
                    index = r.t(delta + n1.index)
                    n3 = AssemblyNode(index, False, str(uuid4()), overhang=True)
                    edata1 = dict(edata)
                    edata2 = dict(edata)
                    edata1["material"] = edata["material"] / 10.0
                    edata2["material"] = edata["material"] / 10.0
                    edata1["span"] = 0

                    gcopy.add_edge(n1, n3, **edata1)
                    gcopy.add_edge(n3, n2, **edata2)
            design.graphs[qk] = gcopy

        result = list(design.optimize().values())[0]
        assembly = result.assemblies[0]
        df = assembly.to_df()
        assert list(df["query_start"]) == [200, 500, 750, 1000, 2000, 3100]
        assert list(df["query_end"]) == [500, 750, 1000, 2000, 3100, 200]
