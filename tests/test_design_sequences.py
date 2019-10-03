import json

from dasi.design.sequence_design import design_primers
from dasi.utils import Region


bases = "AGCT"
bases += bases.lower()

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
