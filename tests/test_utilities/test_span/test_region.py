from dasi.utils import Region


def test_flip_region():
    s1 = "TGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTG"
    r1 = Region(10, 30, len(s1), cyclic=True, direction=1)
    print(r1)

    s2 = s1[::-1]
    r2 = r1.flip()
    print(r2)

    assert r1.get_slice(s1) == r2.get_slice(s2)[::-1]


def test_flip_cyclic_region():
    s1 = "TGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTG"
    r1 = Region(len(s1) - 10, 10, len(s1), cyclic=True, direction=1)
    print(r1)

    s2 = s1[::-1]
    r2 = r1.flip()
    print(r2)

    assert r1.get_slice(s1) == r2.get_slice(s2)[::-1]
