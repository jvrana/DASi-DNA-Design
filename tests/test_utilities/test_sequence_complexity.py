import random

import pytest
from flaky import flaky

from dasi.utils.sequence import DNAStats


def random_seq(length, bases=None):
    if bases is None:
        bases = "AGTC"

    seq = ""
    for _ in range(length):
        i = random.randint(0, len(bases) - 1)
        seq += bases[i]
    return seq


def revcomp(seq):
    r = ""
    d = dict(zip("actgACTG", "tgacTGAC"))
    for s in seq:
        r += d[s]
    return r[::-1]


class TestDNAStats:
    @flaky(max_runs=25, min_passes=25)
    def test_normal_sequence(self):
        seq = random_seq(1000)
        stats = DNAStats(seq, 14, 20, 20)
        print(seq)
        print(stats())
        print(stats.cost())
        assert stats.cost() < 15

        print("OK")
        print(
            DNAStats(
                "CGAGACCACTCGGGACTTCCGGCCATAGCGTACCGTTTTGTGACAAAACCCCCACTCGAACGTGAGAAAACCCTTCTCTTCATGTAATTCCCCCACAGTCCGCGGGTCGGTCAAACCTGGATAAGGTAAAGACTAATATCTAAACCTGCTGGAGAGTCGAACCGCGGTCTTAGGCCCACGCAGAGTGTATGTTATTCGTCTGCCGCTATATCGGTCAACACTAGTTGACGGATAGGAATGTTGGATTAACGCGTCTCCAACGCTGGGATACCCTCGCAAAATTTTCCCGATACTATCCGGAATCTCTAACGCCGTTGGTTTGGGCTCCCAACCACCCGTGAACTTCTAACACGAGAATCACCGCTGGAGCGCGCGCCTTCTCTCAATTTACCTGAGCTTTCGCTTCCTACTTAGCAGAATCGTGAACCTAAATTTTAGCAGCTTCAAGTCAGTTACGCTCGACACTTCCGATTCCAGGTAAAATAACCACTTCTAAGGTTCGTGACTGGTTCTCTATTCAACGCACGCGGTGCCCTCGCGGGTCCTCTGCTGCCGGGAAGCACATGATTGCCAGCTTGTTAAACAACACAAGGTGGCCAATCTCAAACTCGCATAAGCCCTGTTTTTTCTTGCAAGCTGCAACCGAGCATTCCTTCAGTCAGTGGTGGTTTTTCAAAACTATTCCTATGGGTGCTGACACGTGTGTAATTGTTTTCTACTATCTCTCGGTTTATAGCGTAGTTGCCGAGGCTATTGAGTCTCCTTTGCTAATAGCTAAGGTGGAAATTTTTTTTTTTTTGAACCGGGTGAATATACTTGATACATCAATAGCCCCTAGCGTATTGTACCCGTCACGGGCTCAAATACTCTGCCCAGGGCGATACCATGGAAGTTCTCGTAACATACAATGGATCTGGGCCGTCATCGCTTGATGCTCTAGAAGAAAAAGCAGAGACCGGCCATTACCGCGTCAACTAACACGCCTCAGGCCGGGGTTAACACTAGGTGTGT",
                14,
                20,
                20,
            )()
        )

    @flaky(max_runs=25, min_passes=25)
    @pytest.mark.parametrize("kmer", [(25, 20), (25, 24), (25, 25)])
    def test_repeats(self, kmer):
        repeat = random_seq(kmer[0])
        seq = random_seq(1000) + repeat + random_seq(500) + repeat + random_seq(1000)

        stats = DNAStats(seq, kmer[1], 20, 20)
        print(stats())
        assert stats()["n_repeats"] > 0
        assert stats()["n_hairpins"] == 0

    @flaky(max_runs=25, min_passes=25)
    @pytest.mark.parametrize("kmer", [(25, 20), (25, 24), (25, 25)])
    def test_hairpins(self, kmer):
        hairpin = random_seq(kmer[0])
        seq = (
            random_seq(1000)
            + hairpin
            + random_seq(500)
            + revcomp(hairpin)
            + random_seq(1000)
        )

        stats = DNAStats(seq, 14, 20, kmer[1])
        print(stats())
        assert stats()["n_repeats"] == 0
        assert stats()["n_hairpins"] > 0

        print(stats.cost(1000, 1500))
        print(stats.cost(None, None))

    @flaky(max_runs=3, min_passes=3)
    def test_partitioner(self):
        repeat = random_seq(40)
        seq = random_seq(2000) + repeat + random_seq(100) + revcomp(repeat)

        stats = DNAStats(seq, 14, 20, 20)

        partitions = stats.partition(10, overlap=25, stopping_threshold=10)
        p = partitions[0]
        print(p)

        safe_start = 2000 + 40
        safe_end = 2000 + 40 + 100
        assert p["index_1"][1] > safe_start
        assert p["index_1"][1] < safe_end

        assert p["index_2"][0] >= safe_start - 14
        assert p["index_2"][0] < safe_end
        assert p["cost"] < 10
        assert p["cost"] < stats.cost()

    @flaky(max_runs=3, min_passes=3)
    def test_partitioner_cached(self):
        repeat = random_seq(40)
        seq = random_seq(2000) + repeat + random_seq(100) + revcomp(repeat)

        stats = DNAStats(seq, 13, 20, 20)

        partitions = stats.partition(10, overlap=25, stopping_threshold=10)
        partitions = stats.partition(
            10, overlap=25, stopping_threshold=10, i=10, j=None
        )
        p = partitions[0]
        print(p)

        safe_start = 2000 + 40
        safe_end = 2000 + 40 + 100
        assert p["index_1"][1] > safe_start
        assert p["index_1"][1] < safe_end

        assert p["index_2"][0] >= safe_start - 14
        assert p["index_2"][0] < safe_end
        assert p["cost"] < 10
        assert p["cost"] < stats.cost()
