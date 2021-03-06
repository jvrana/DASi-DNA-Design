import functools
import itertools
import operator
import random
from typing import List
from typing import Union

import pytest
from flaky import flaky

from dasi.utils.sequence.sequence_complexity import count_misprimings_in_amplicon
from dasi.utils.sequence.sequence_complexity import DNAStats


##########################################
# Utility methods
##########################################


def random_seq(length, bases=None):
    """Produce a randomized sequence."""
    if bases is None:
        bases = "AGTC"

    seq = ""
    for _ in range(length):
        i = random.randint(0, len(bases) - 1)
        seq += bases[i]
    return seq


def revcomp(seq):
    """Reverse complement a sequence."""
    r = ""
    d = dict(zip("actgACTG", "tgacTGAC"))
    for s in seq:
        r += d[s]
    return r[::-1]


def cross_product(*args: List[Union[List[str], str]], operation=operator.add):
    """Develop sequences by taking a cross product of a list of lists."""
    list_of_lists = []
    for arg in args:
        if isinstance(arg, list):
            list_of_lists.append(arg)
        else:
            list_of_lists.append([arg])

    sequences = []
    for x in itertools.product(*list_of_lists):
        if operation:
            sequences.append(functools.reduce(operation, x))
    return sequences


def test_cross_product():
    seqs = cross_product(["1"], ["1", "2"], ["1", "2", "3"], "1")
    assert len(seqs) == 6
    for s in seqs:
        print(s)


def replace(to_replace, replace_with, a, b, c, rc):
    """Replace a section of the sequence `to_replace` with a secion of the
    sequence `replace_with`

    :param to_replace: position to replace sequence
    :param replace_with: indices of the sequence to replace with replace_with[a:b]
    :param a: starting index of the replacement sequence
    :param b: ending index of the replacement sequence
    :param c: index to replace sequence.
    :param rc: if True, replaces with reverse_complement
    :return:
    """
    x1 = replace_with[a:b]
    if rc:
        x1 = revcomp(x1)
    return to_replace[:c] + x1 + to_replace[c + len(x1) :]


def replace_rc_pairs(to_replace, replace_with, a, b, c):
    """Return `replace(..., rc=False)` and `replace(..., rc=True)`"""
    x1 = replace_with[a:b]
    x2 = revcomp(x1)
    s1 = to_replace[:c] + x1 + to_replace[c + len(x1) :]
    s2 = to_replace[:c] + x2 + to_replace[c + len(x2) :]
    assert len(s1) == len(s2)
    return s1, s2


def gen_seq(
    i, j, length, window, insert_loc, repeat_i, repeat_l, rc, to_replace, replace_with
):
    repeat_j = repeat_l + repeat_i

    left_flank = random_seq(i)
    right_flank = random_seq(length - j)
    amplicon_left_window = random_seq(window)
    amplicon_internal = random_seq(j - i - window * 2)
    amplicon_right_window = random_seq(window)

    seqs = {
        "left_window": amplicon_left_window,
        "internal": amplicon_internal,
        "right_window": amplicon_right_window,
        "left_flank": left_flank,
        "right_flank": right_flank,
    }

    # add the repeat
    seqs[to_replace] = replace(
        seqs[to_replace], seqs[replace_with], repeat_i, repeat_j, insert_loc, rc
    )

    amplicon = seqs["left_window"] + seqs["internal"] + seqs["right_window"]

    seq = seqs["left_flank"] + amplicon + seqs["right_flank"]
    assert len(seq) == length, "testing sequence was constructed incorrectly"
    return seq


def gen_seq_with_rc_pairs(
    i, j, length, window, insert_loc, repeat_i, repeat_l, to_replace, replace_with
):
    repeat_j = repeat_l + repeat_i

    left_flank = "N" * i  # random_seq(i)
    right_flank = "N" * (length - j)  # random_seq(l - j)
    amplicon_left_window = random_seq(window)
    amplicon_internal = "N" * (j - i - window * 2)  # random_seq(j - i - window * 2)
    amplicon_right_window = random_seq(window)

    seqs = {
        "left_window": amplicon_left_window,
        "internal": amplicon_internal,
        "right_window": amplicon_right_window,
        "left_flank": left_flank,
        "right_flank": right_flank,
    }

    seqs[to_replace] = list(
        replace_rc_pairs(
            seqs[to_replace], seqs[replace_with], repeat_i, repeat_j, insert_loc
        )
    )

    sequences = cross_product(
        seqs["left_flank"],
        seqs["left_window"],
        seqs["internal"],
        seqs["right_window"],
        seqs["right_flank"],
    )

    for s in sequences:
        if not len(s) == length:
            raise ValueError("{} != {}".format(len(s), length))

    return sequences


# DNAStats tests
class TestDNAStats:
    """Test the DNAStats instance."""

    def test_case_insensitive(self):
        seq = random_seq(1000)
        stats1 = DNAStats(seq.lower(), 20, 20, 20)
        stats2 = DNAStats(seq.upper(), 20, 20, 20)
        print(stats1.cost())
        assert stats1.cost() == stats2.cost()

    def test_fwd_signatures(self):
        seq = "N" * 100 + random_seq(100) + "N" * 120
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        assert set(stats.fwd_signatures[:70]) == {0.0}
        assert stats.fwd_signatures[71] != 0.0
        assert stats.fwd_signatures[72] != 0.0

    def test_rev_signatures(self):
        seq = "N" * 100 + random_seq(100) + "N" * 101
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        assert set(stats.rev_signatures[:70]) == {0.0}
        assert stats.rev_signatures[71] != 0.0
        assert stats.rev_signatures[72] != 0.0

    @flaky(max_runs=25, min_passes=25)
    def test_normal_sequence(self):
        seq = random_seq(1000)
        stats = DNAStats(seq, 14, 20, 20)
        print(seq)
        print(stats())
        print(stats.cost())
        assert stats.cost() < 15

        print(
            DNAStats(
                "CGAGACCACTCGGGACTTCCGGCCATAGCGTACCGTTTTGTGACAAA"
                "ACCCCCACTCGAACGTGAGAAAACCCTTCTCTTCATGTAATTCCCCCACA"
                "GTCCGCGGGTCGGTCAAACCTGGATAAGGTAAAGACTAATATCTAAACCT"
                "GCTGGAGAGTCGAACCGCGGTCTTAGGCCCACGCAGAGTGTATGTTA"
                "TTCGTCTGCCGCTATATCGGTCAACACTAGTTGACGGATAGGAATGTTGG"
                "ATTAACGCGTCTCCAACGCTGGGATACCCTCGCAAAATTTTCCCGAT"
                "ACTATCCGGAATCTCTAACGCCGTTGGTTTGGGCTCCCAACCACCCGTG"
                "AACTTCTAACACGAGAATCACCGCTGGAGCGCGCGCCTTCTCTCAATT"
                "TACCTGAGCTTTCGCTTCCTACTTAGCAGAATCGTGAACCTAAATTTTA"
                "GCAGCTTCAAGTCAGTTACGCTCGACACTTCCGATTCCAGGTAAAATA"
                "ACCACTTCTAAGGTTCGTGACTGGTTCTCTATTCAACGCACGCGGTGCCC"
                "TCGCGGGTCCTCTGCTGCCGGGAAGCACATGATTGCCAGCTTGTTAA"
                "ACAACACAAGGTGGCCAATCTCAAACTCGCATAAGCCCTGTTTTTTCTTG"
                "CAAGCTGCAACCGAGCATTCCTTCAGTCAGTGGTGGTTTTTCAAAAC"
                "TATTCCTATGGGTGCTGACACGTGTGTAATTGTTTTCTACTATCTCTCG"
                "GTTTATAGCGTAGTTGCCGAGGCTATTGAGTCTCCTTTGCTAATAGCT"
                "AAGGTGGAAATTTTTTTTTTTTTGAACCGGGTGAATATACTTGATACAT"
                "CAATAGCCCCTAGCGTATTGTACCCGTCACGGGCTCAAATACTCTGCC"
                "CAGGGCGATACCATGGAAGTTCTCGTAACATACAATGGATCTGGGCCGT"
                "CATCGCTTGATGCTCTAGAAGAAAAAGCAGAGACCGGCCATTACCGCG"
                "TCAACTAACACGCCTCAGGCCGGGGTTAACACTAGGTGTGT",
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
    @pytest.mark.parametrize("length__kmer", [(25, 20), (25, 24), (25, 25)])
    def test_hairpins(self, length__kmer):
        hairpin_length, kmer = length__kmer
        hairpin = random_seq(hairpin_length)

        # make a sequence with a hairpin
        seq = (
            random_seq(1000)
            + hairpin
            + random_seq(500)
            + revcomp(hairpin)
            + random_seq(1000)
        )

        # look for hairpins of size kmer
        stats = DNAStats(seq, 20, 20, kmer)

        assert stats()["n_repeats"] == 0
        assert stats()["n_hairpins"] > 0

        print(stats.cost(1000, 1500))
        print(stats.cost(None, None))

    def test_count_misprimings_from_slice_case1(self):
        repeat = random_seq(30)
        seq = "T" * 100 + repeat + "T" * 100
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        n = stats.count_repeats_from_slice(100, 130)
        assert n == 0

    def test_count_misprimings_from_slice_case2(self):
        repeat = random_seq(30)
        seq = "N" * 100 + repeat + "N" * 100 + repeat + "N" * 107
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        n = stats.count_repeats_from_slice(100, 130)
        assert n == 1

    def test_count_misprimings_from_slice_case3(self):
        """repeats are very near the edge."""
        repeat = random_seq(30)
        seq = "N" + repeat + "N" * 100 + revcomp(repeat) + "N"
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        n = stats.count_repeats_from_slice(1, 31)
        assert n == 1

    @pytest.mark.parametrize("ij", [(100, 130), (160, 190), (260, 290)])
    def test_count_misprimings_from_slice_case4(self, ij):
        """Here we have a sequence of Ns with a predicted hairpin at indices.

        [100:130] and [160:190] and [260:290].

        Evaluating at any of these indices should return exactly 2
        sequences.
        """
        repeat = random_seq(30)
        i, j = ij
        seq = (
            "N" * 100
            + repeat
            + "N" * 30
            + repeat
            + "N" * 100
            + revcomp(repeat)
            + "N" * 107
        )
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        n = stats.count_repeats_from_slice(i, j)
        assert n == 2

    def test_count_misprimings_from_slice_case5(self):
        """Here we have a sequence of Ns with a predicted hairpin at indices.

        [100:130] and [200:230].
        """
        repeat = random_seq(30)
        seq = "N" * 100 + repeat + "N" * 100 + revcomp(repeat) + "N" * 107
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        n = stats.count_repeats_from_slice(100, 130)
        assert n == 1

    @pytest.mark.parametrize("dist_from_left", [0, 1])
    @pytest.mark.parametrize("dist_from_right", [0, 1])
    @pytest.mark.parametrize("left_or_right", ["left", "right"])
    @pytest.mark.parametrize("rc", [True, False])
    def test_count_misprimings_from_slice_edge_cases(
        self, dist_from_left, dist_from_right, left_or_right, rc
    ):
        """Here one of the repeats is on the 3' (right) end of the sequence."""
        repeat = random_seq(30)
        r1 = repeat
        r2 = repeat

        if left_or_right == "left" and rc:
            r1 = revcomp(r1)
        elif left_or_right == "right" and rc:
            r2 = revcomp(r2)

        seq = "N" * dist_from_left + r1 + "N" * 100 + r2 + "N" * dist_from_right
        stats = DNAStats(seq, 1, 1, hairpin_window=30)
        print(seq)
        print(seq[dist_from_left:dist_from_left+30])
        if left_or_right == "left":
            n = stats.count_repeats_from_slice(dist_from_left, dist_from_left + 30)
            assert n == 1
        else:
            n = stats.count_repeats_from_slice(
                dist_from_left + 30 + 100, dist_from_left + 30 + 100 + 30
            )
            assert n == 1


class TestPositiveRCExamples:
    """Tests that reverse_complemented repeats results in the name number of
    misprimings."""

    @flaky(max_runs=10, min_passes=10)
    @pytest.mark.parametrize(
        "to_replace",
        ["left_flank", "right_flank", "internal"],
        ids=["left_flank", "right_flank", "internal"],
    )
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["left_window", "right_window"],
    )
    @pytest.mark.parametrize("cyclic", [True, False], ids=["cyclic", "linear"])
    @pytest.mark.parametrize("repeat_l", [20, 21, 22])
    def test_cyclic_and_linear_same_misprimings(
        self, to_replace, replace_with, repeat_l, cyclic
    ):
        """Should not matter if the sequence being replace is reverse-
        complement or not."""

        window = 30
        i, j = 100, 250
        kwargs = dict(
            i=i,
            j=j,
            length=500,
            window=window,
            insert_loc=10,
            repeat_i=5,
            repeat_l=repeat_l,
            to_replace=to_replace,
            replace_with=replace_with,
        )
        seq1, seq2 = gen_seq_with_rc_pairs(**kwargs)
        min_primer_anneal = 20
        n1 = count_misprimings_in_amplicon(
            seq1,
            i,
            j,
            min_primer_anneal=min_primer_anneal,
            max_primer_anneal=window,
            cyclic=cyclic,
        )
        n2 = count_misprimings_in_amplicon(
            seq2,
            i,
            j,
            min_primer_anneal=min_primer_anneal,
            max_primer_anneal=window,
            cyclic=cyclic,
        )
        print(n1, n2)
        assert n1 == n2
        assert n1 == repeat_l - min_primer_anneal + 1


@pytest.mark.parametrize("ij", [(100, 300), (500, 800), (0, 500), (800, 1000)])
def test_count_misprimings_in_amplicon(ij):
    seq = random_seq(1000)
    count_misprimings_in_amplicon(
        seq, ij[0], ij[1], min_primer_anneal=12, max_primer_anneal=30
    )


@pytest.mark.parametrize(
    "repeat_length__kmer__num_misprimings",
    [(23, 20, 4), (22, 20, 3), (20, 20, 1), (12, 20, 0)],
)
@pytest.mark.parametrize("rc", [False, True], ids=["", "reverse_complement"])
@pytest.mark.parametrize("repeat_i", [0])
class TestPositiveExamples:
    """
    Tests if we can find the repeats. All tests are parametrized in the following
    way:

    - n_misprimings - the minimum number of misprimings we expect.
    - kmer - the size of the kmer we are looking for in the algorith
    - repeat_length - the size of the repeat region we are placing into the testing
            sequences. If this size is less than the size of the kmer, we should not
            be able to detect it.
    """

    N_RUNS = 3
    LENGTH = 1000
    WINDOW = 30

    @flaky(max_runs=N_RUNS, min_passes=N_RUNS)
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["with_left_window", "with_right_window"],
    )
    @pytest.mark.parametrize("to_replace", ["internal"], ids=["replace_internal"])
    @pytest.mark.parametrize("ij", [(200, 800), (150, 850)])
    @pytest.mark.parametrize("cyclic", [False, True], ids=["linear", "cyclic"])
    def test_find_repeat_inside_internal_amplicon(
        self,
        repeat_length__kmer__num_misprimings,
        rc,
        to_replace,
        replace_with,
        repeat_i,
        ij,
        cyclic,
    ):
        """Tests if we can find the repeat from within the pcr product from the
        left window."""
        repeat_length, kmer, n_misprimings = repeat_length__kmer__num_misprimings
        i, j = ij
        window = self.WINDOW
        length = self.LENGTH
        seq = gen_seq(
            i=i,
            j=j,
            length=length,
            window=window,
            insert_loc=100,
            repeat_i=repeat_i,
            repeat_l=repeat_length,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )
        n = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window, cyclic=cyclic
        )
        assert n >= n_misprimings
        if n_misprimings == 0:
            assert n == 0

    @flaky(max_runs=N_RUNS, min_passes=N_RUNS)
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["with_left_window", "with_right_window"],
    )
    @pytest.mark.parametrize(
        "to_replace",
        ["left_flank", "right_flank"],
        ids=["replace_left_flank", "replace_right_flank"],
    )
    @pytest.mark.parametrize("ij", [(200, 800), (150, 850)])
    @pytest.mark.parametrize("cyclic", [False, True], ids=["linear", "cyclic"])
    def test_find_repeat_inside_external_to_amplicon(
        self,
        repeat_length__kmer__num_misprimings,
        rc,
        to_replace,
        replace_with,
        repeat_i,
        ij,
        cyclic,
    ):
        """Tests if we can find the repeat from within the pcr product from the
        left window."""
        repeat_length, kmer, n_misprimings = repeat_length__kmer__num_misprimings
        i, j = ij
        window = self.WINDOW
        length = self.LENGTH
        seq = gen_seq(
            i=i,
            j=j,
            length=length,
            window=window,
            insert_loc=100,
            repeat_i=repeat_i,
            repeat_l=repeat_length,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )
        n = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window, cyclic=cyclic
        )
        assert n >= n_misprimings
        if n_misprimings == 0:
            assert n == 0

    @flaky(max_runs=N_RUNS, min_passes=N_RUNS)
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["with_left_window", "with_right_window"],
    )
    @pytest.mark.parametrize("to_replace", ["internal"], ids=["replace_internal"])
    @pytest.mark.parametrize("ij", [(200, 800), (150, 850)])
    def test_linear_and_cyclic_same_bindings(
        self,
        repeat_length__kmer__num_misprimings,
        rc,
        to_replace,
        replace_with,
        repeat_i,
        ij,
    ):
        """Tests if n bindings is the same between cyclic and linear
        templates."""
        repeat_length, kmer, n_misprimings = repeat_length__kmer__num_misprimings
        i, j = ij
        window = self.WINDOW
        length = self.LENGTH
        seq = gen_seq(
            i=i,
            j=j,
            length=length,
            window=window,
            insert_loc=100,
            repeat_i=repeat_i,
            repeat_l=repeat_length,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )
        n1 = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window, cyclic=False
        )
        n2 = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window, cyclic=True
        )
        assert n1 == n2


class TestPositiveCyclicExamples:
    @flaky(max_runs=10, min_passes=10)
    @pytest.mark.parametrize(
        "to_replace",
        ["left_flank", "right_flank", "internal"],
        ids=["left_flank", "right_flank", "internal"],
    )
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["left_window", "right_window"],
    )
    @pytest.mark.parametrize("rc", [True, False], ids=["rc", ""])
    def test_cyclic_and_linear_same_misprimings(self, to_replace, replace_with, rc):
        """Should not matter if the sequence is reindex."""

        window = 30
        i, j = 200, 800
        seq = gen_seq(
            i=i,
            j=j,
            length=1000,
            window=window,
            insert_loc=100,
            repeat_i=0,
            repeat_l=25,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )

        n1 = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=20, max_primer_anneal=window, cyclic=False
        )
        n2 = count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=20, max_primer_anneal=window, cyclic=True
        )
        print(n1, n2)
        assert n1 == n2

    @pytest.mark.parametrize(
        "to_replace",
        ["left_flank", "right_flank", "internal"],
        ids=["left_flank", "right_flank", "internal"],
    )
    @pytest.mark.parametrize(
        "replace_with",
        ["left_window", "right_window"],
        ids=["left_window", "right_window"],
    )
    def test_cyclic_indices(self, to_replace, replace_with):
        length = 1000
        i = 200
        j = 800
        window = 30

        seq1, seq2 = gen_seq_with_rc_pairs(
            i=i,
            j=j,
            length=length,
            window=window,
            insert_loc=10,
            repeat_i=0,
            repeat_l=20,
            to_replace=to_replace,
            replace_with=replace_with,
        )

        # linear misprimings
        n1 = count_misprimings_in_amplicon(
            seq1, i, j, min_primer_anneal=20, max_primer_anneal=window, cyclic=False
        )

        # reindex the plasmid 100bp in the i:j slice
        shift = 100
        delta = i + shift
        seq3 = seq1[delta:] + seq1[:delta]
        i = length - shift
        j = j - delta

        # cyclic misprimings
        n2 = count_misprimings_in_amplicon(
            seq3, i, j, min_primer_anneal=20, max_primer_anneal=window, cyclic=True
        )

        assert n1 == n2
        assert n1 == 1

    def test_cyclic_indices_raises_index_error(self):
        kwargs = dict(
            seq=random_seq(2000),
            i=900,
            j=200,
            min_primer_anneal=20,
            max_primer_anneal=30,
        )
        with pytest.raises(IndexError):
            count_misprimings_in_amplicon(**kwargs, cyclic=False)
        count_misprimings_in_amplicon(**kwargs, cyclic=True)


class TestNegativeExamples:
    """Tests situations in which we will NOT count repeats."""

    @pytest.mark.parametrize("to_replace", ["left_flank", "right_flank", "internal"])
    @pytest.mark.parametrize("replace_with", ["left_window", "right_window"])
    @pytest.mark.parametrize("rc", [True, False])
    @pytest.mark.parametrize("window_results", [(15, False), (30, True)])
    def test_window(self, to_replace, replace_with, rc, window_results):
        """We expect that if the window is too small, we will not count any
        repeats."""
        i, j = 150, 850
        repeat_i = 0
        length = 1000
        kmer = 20
        window, expect = window_results
        seq = gen_seq(
            i=i,
            j=j,
            length=length,
            window=window,
            insert_loc=100,
            repeat_i=repeat_i,
            repeat_l=20,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )
        if not expect:
            count_misprimings_in_amplicon(
                seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window
            ) == 0
        else:
            count_misprimings_in_amplicon(
                seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window
            ) > 0

        count_misprimings_in_amplicon(
            seq, i, j, min_primer_anneal=20, max_primer_anneal=window, cyclic=False
        )

    @pytest.mark.parametrize(
        "replacements",
        [
            ("left_flank", "right_flank"),
            ("right_flank", "left_flank"),
            ("internal", "right_flank"),
            ("internal", "left_flank"),
            ("left_flank", "internal"),
            ("right_flank", "internal"),
        ],
    )
    @pytest.mark.parametrize("rc", [True, False])
    @pytest.mark.parametrize("kmer", [18, 20, 22])
    @pytest.mark.parametrize("repeat_length", [18, 20, 22])
    @pytest.mark.parametrize("repeat_i", [0, 10])
    def test_external_repeat_has_no_misprimings(
        self, rc, replacements, repeat_length, kmer, repeat_i
    ):
        """If a sequence has a repeat not in the expected primer window, then
        we will not count it as a repeat."""
        to_replace, replace_with = replacements
        n_misprimings = 0
        i = 200
        j = 800
        window = 30
        seq = gen_seq(
            i=i,
            j=j,
            length=1000,
            window=window,
            insert_loc=100,
            repeat_i=repeat_i,
            repeat_l=repeat_length,
            to_replace=to_replace,
            replace_with=replace_with,
            rc=rc,
        )
        assert (
            count_misprimings_in_amplicon(
                seq, i, j, min_primer_anneal=kmer, max_primer_anneal=window
            )
            >= n_misprimings
        )


def test_hash1():
    s1 = random_seq(1000)
    f = functools.partial(
        DNAStats, repeat_window=14, stats_window=20, hairpin_window=20
    )
    stats1 = f(s1)
    stats2 = f(s1)
    stats3 = f(random_seq(1000))
    print(hash(stats1))
    print(hash(stats2))
    print(hash(stats3))
    assert hash(stats1) == hash(stats2)
    assert not hash(stats1) == hash(stats3)


@pytest.mark.parametrize("key", ["repeat_window", "stats_window", "hairpin_window"])
def test_hash2(key):
    s1 = random_seq(1000)
    kwargs = {"repeat_window": 20, "stats_window": 20, "hairpin_window": 20}
    kwargs2 = dict(kwargs)
    kwargs2[key] += 1
    stats1 = DNAStats(s1, **kwargs)
    stats2 = DNAStats(s1, **kwargs)
    stats3 = DNAStats(s1, **kwargs2)
    assert hash(stats1) == hash(stats2)
    assert not hash(stats1) == hash(stats3)


@pytest.mark.parametrize(
    "key", ["gc_content_threshold", "at_content_threshold", "base_percentage_threshold"]
)
def test_hash3(key):
    s1 = random_seq(1000)
    kwargs = {
        "repeat_window": 20,
        "stats_window": 20,
        "hairpin_window": 20,
        "gc_content_threshold": 0.8,
        "at_content_threshold": 0.8,
        "base_percentage_threshold": 0.8,
    }
    kwargs2 = dict(kwargs)
    kwargs2[key] += 0.1
    stats1 = DNAStats(s1, **kwargs)
    stats2 = DNAStats(s1, **kwargs)
    stats3 = DNAStats(s1, **kwargs2)
    assert hash(stats1) == hash(stats2)
    assert not hash(stats1) == hash(stats3)


# TODO: actually test that underlying data is the same for the view
def test_view():
    seq = random_seq(1000)
    stats = DNAStats(seq, 14, 20, 20)
    stats2 = stats.view(slice(None, None))
    assert stats is not stats2
    assert stats.cost(1, 1000) == stats.cost(1, 1000)
    print(stats)
    print(stats2)


def test_copy():
    seq = random_seq(1000)
    stats = DNAStats(seq, 14, 20, 20)
    stats2 = stats.copy(slice(None, None))
    assert stats is not stats2
    assert stats.cost(1, 1000) == stats.cost(1, 1000)
    print(stats)
    print(stats2)


@pytest.mark.parametrize("index", [30, 400, None])
def test_slice(index):
    seq = random_seq(1000)
    stats = DNAStats(seq, 14, 20, 20)
    stats2 = stats[:index]
    if index is None:
        assert len(stats2) == len(stats)
    else:
        assert len(stats2) == index
    print(stats.cost())
    print(stats2.cost())
