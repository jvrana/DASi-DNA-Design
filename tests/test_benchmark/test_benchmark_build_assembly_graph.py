from shoestring import AssemblyGraphBuilder, AlignmentContainer, Alignment, Constants, Region
import pytest
from shoestring.utils import make_async

fake_container = AlignmentContainer({})

builder = AssemblyGraphBuilder(fake_container)
builder.build_assembly_graph()



def make_fake_alignment(a, b):
    alignment = Alignment(
        Region(a, b, 10000),
        Region(a, b, 10000),
        type=Constants.PCR_PRODUCT,
        query_key="",
        subject_key="",
    )
    fake_container.alignments.append(alignment)
    return alignment


@pytest.mark.parametrize('r', [
    (0, 100, 30, 110),
    (0, 100, 75, 175),
    (0, 100, 20, 80),
    (0, 100, 200, 300)
], ids=[
    "overlap left",
    "overlap right",
    "contains",
    "no overlap"
])
def test_add_edge(benchmark, r):
    make_fake_alignment(*r[:2])
    make_fake_alignment(*r[2:])

    builder = AssemblyGraphBuilder(fake_container)
    groups = fake_container.alignment_groups
    builder.build_assembly_graph()

    benchmark(builder.add_edge, groups[0], groups[1])


@pytest.mark.parametrize('num_threads', [1, 10])
def test_add_multiple_edges(benchmark, num_threads):
    r = (0, 100, 30, 110)
    make_fake_alignment(*r[:2])
    make_fake_alignment(*r[2:])

    groups = fake_container.alignment_groups


    @make_async(num_threads)
    def add_edges(edges):
        for g1, g2 in edges:
            builder.add_edge(g1, g2)

    edges = [(groups[0], groups[1])] * 10000

    benchmark(add_edges, edges)