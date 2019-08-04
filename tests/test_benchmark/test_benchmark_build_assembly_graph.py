from shoestring import AssemblyGraphBuilder, AlignmentContainer, Alignment, Constants, Region
import pytest


fake_container = AlignmentContainer({})


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
