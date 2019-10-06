import json
import os
from typing import Dict

import dill
import pytest
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.alignments import AlignmentGroup
from dasi.cost import SpanCost
from dasi.design import Design
from dasi.design.sequence_design import design_primers
from dasi.utils import Region

gfp = "ATGGTCTCTAAGGGTGAAGAATTGTTCACCGGTGTCGTCCCAATCTTGGTCGAATTGGACGGGGACGTCAACGGTCACAAGTTCTCTGTCTCTGGTGAAGGTGAAGGTGACGCTACCTACGGTAAGTTGACCTTGAAGTTCATCTGTACCACCGGTAAGTTGCCAGTCCCATGGCCAACCTTGGTCACCACCTTCGGTTACGGTGTCCAATGTTTCGCTAGATACCCAGACCACATGAAGCAACACGACTTCTTCAAGTCTGCTATGCCAGAAGGTTACGTCCAAGAAAGAACCATCTTCTTCAAGGACGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTCGAAGGTGACACCTTGGTCAACAGAATCGAATTGAAGGGTATCGACTTCAAGGAAGACGGTAACATCTTGGGTCACAAGTTGGAATACAACTACAACTCTCACAACGTCTACATCATGGCTGACAAGCAAAAGAACGGTATCAAGGTCAACTTCAAGATCAGACACAACATCGAAGACGGTTCTGTCCAATTGGCTGACCACTACCAACAAAACACCCCAATCGGTGACGGTCCAGTCTTGTTGCCAGACAACCACTACTTGTCTACCCAATCTGCTTTGTCTAAGGACCCAAACGAAAAGAGAGACCACATGGTCTTGTTGGAATTCGTCACCGCTGCTGGTATCACCCACGGTATGGACGAATTGTACAAGTAA"


def test_region_invert():
    """Check Region.invert.

    This is essential to the primer design algorithm.
    """
    r = Region(10, 20, 50, cyclic=True)
    r1, r2 = r.invert()
    assert r2 is None
    assert list(r1) == list(range(20, 50)) + list(range(10))


def test_primer_design():
    region = Region(100, 300, len(gfp), cyclic=True)
    pairs, explain = design_primers(gfp, region, None, None)
    print(json.dumps(pairs, indent=1))

    for pair in pairs.values():
        assert pair["LEFT"]["location"][0] == 100
        assert pair["RIGHT"]["location"][0] == 299


def test_primer_design2():
    i = len(gfp) - 50
    j = 100
    region = Region(i, j, len(gfp), cyclic=True)
    pairs, explain = design_primers(gfp, region, None, None)
    # print(json.dumps(pairs, indent=1))
    for pair in pairs.values():
        print(pair["LEFT"]["location"])
        assert pair["LEFT"]["location"][0] == i
        assert pair["RIGHT"]["location"][0] == j - 1


def pkl_results(here, paths, query, span_cost):
    path = "results.pkl"

    if os.path.isfile(path):
        with open(path, "rb") as f:
            return dill.load(f)
    else:
        primers = make_linear(load_fasta_glob(paths["primers"]))
        templates = load_genbank_glob(paths["templates"]) + load_genbank_glob(
            paths["registry"]
        )

        query_path = os.path.join(here, "data/test_data/genbank/designs", query)
        queries = make_circular(load_genbank_glob(query_path))

        design = Design(span_cost=span_cost)

        design.add_materials(primers=primers, templates=templates, queries=queries)

        design.compile()

        assert len(design.graphs) == len(queries)
        assert len(design.graphs) == 1

        results = design.optimize()
        with open(path, "wb") as f:
            dill.dump(results, f)
        return results


@pytest.mark.parametrize("query", ["pmodkan-ho-pact1-z4-er-vpr.gb"])
def test_design_with_primers(here, paths, query, span_cost):
    query_path = "pmodkan-ho-pact1-z4-er-vpr.gb"
    results = pkl_results(here, paths, query_path, span_cost)

    # TEST HERE
    from dasi.design.sequence_design import design_primers
    from dasi.utils import NumpyDataFrame
    from dasi.design import Assembly, AssemblyNode
    import numpy as np

    def edata_to_npdf(edata: dict, span_cost: SpanCost) -> NumpyDataFrame:
        npdf = span_cost.cost(np.array([edata["span"]]), edata["type_def"].design)
        return npdf[0]

    def no_none_or_nan(*i):
        for _i in i:
            if _i is not None and not np.isnan(_i):
                return _i

    def get_primer_extensions(graph, n1, n2, cyclic=True):
        successors = list(graph.successors(n2))
        if successors:
            sedge = graph[n2][successors[0]]
            right = edata_to_npdf(sedge, span_cost)
            r1 = right.data["rprimer_right_ext"]
            r2 = right.data["right_ext"]
            right_ext = no_none_or_nan(r2, r1)
        elif cyclic:
            raise Exception
        else:
            right_ext = 0

        predecessors = list(graph.predecessors(n1))
        if predecessors:
            pedge = graph[predecessors[0]][n1]
            left = edata_to_npdf(pedge, span_cost)
            l1 = left.data["lprimer_left_ext"]
            l2 = left.data["left_ext"]
            left_ext = no_none_or_nan(l2, l1)
        elif cyclic:
            raise Exception
        else:
            left_ext = 0
        return left_ext, right_ext

    def design_pcr_product_primers(
        assembly: Assembly,
        n1: AssemblyNode,
        n2: AssemblyNode,
        span_cost: SpanCost,
        seqdb,
    ):
        graph = assembly.graph

        edge = graph[n1][n2]
        moltype = edge["type_def"]

        if moltype.use_direct:
            print("USE DIRECTLY")
        elif moltype.synthesize:
            gene = edata_to_npdf(edge, span_cost)
            gene_size = gene.data["gene_size"]
            lshift = gene.data["lshift"]
            assert not np.isnan(gene_size)
            assert not np.isnan(lshift)
        else:
            left_ext, right_ext = get_primer_extensions(graph, n1, n2)

            # contains information about templates and queries
            alignment_groups = edge["groups"]
            group = alignment_groups[0]

            lkey, rkey = None, None
            if moltype.design == (1, 0):
                tkey, rkey = group.subject_keys
                region = group.alignments[0].subject_region
            elif moltype.design == (0, 1):
                lkey, tkey = group.subject_keys
            elif moltype.design == (1, 1):
                assert isinstance(group, AlignmentGroup)
                tkey = group.subject_keys[0]
                region = group.alignments[0].subject_region
            else:
                raise Exception("Edge type not understood.")

            trecord = seqdb[tkey]
            tseq = str(trecord.seq)
            if rkey:
                rrecord = seqdb[rkey]
                rseq = str(rrecord.seq)
            else:
                rrecord = None
                rseq = None

            if lkey:
                lrecord = seqdb[lkey]
                lseq = str(lrecord.seq)
            else:
                lrecord = None
                lseq = None
            print("DESIGNING PRIMERS!")
            pairs, explain = design_primers(tseq, region, lseq, rseq)
            print(explain)
            assert pairs

    result = results[list(results)[0]]
    assembly = result.assemblies[0]
    seqdb = result.container.seqdb
    for n1, n2, edata in assembly.edges():
        design_pcr_product_primers(assembly, n1, n2, span_cost, seqdb)

    d = span_cost(1000, (0, 0))
    print(d.to_df().T)
