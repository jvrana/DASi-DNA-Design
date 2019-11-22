import random

from flaky import flaky

from dasi.utils.sequence_complexity import complexity
from dasi.utils.sequence_complexity import complexity_score


def random_seq(length, bases=None):
    if bases is None:
        bases = "AGTC"

    seq = ""
    for _ in range(length):
        i = random.randint(0, len(bases) - 1)
        seq += bases[i]
    return seq


@flaky(max_runs=10, min_passes=9)
def test_complexity_score():
    seq = random_seq(1000)
    print(complexity(seq))
    assert complexity_score(seq) < 10


@flaky(max_runs=10, min_passes=9)
def test_gc_complexity():
    seq1 = random_seq(1000)
    seq2 = random_seq(1000, bases="GGCCAT")
    assert complexity_score(seq2) > complexity_score(seq1)


@flaky(max_runs=10, min_passes=9)
def test_at_polymetrix_streak():
    seq1 = random_seq(1000)
    length = 15
    seq2 = seq1[:500] + random_seq(length, bases="AT") + seq1[500 - length :]
    assert complexity_score(seq2) > complexity_score(seq1)


@flaky(max_runs=10, min_passes=9)
def test_gc_polymetrix_streak():
    seq1 = random_seq(1000)
    length = 10
    seq2 = seq1[:500] + random_seq(length, bases="GC") + seq1[500 - length :]
    assert complexity_score(seq2) > complexity_score(seq1)


@flaky(max_runs=10, min_passes=9)
def test_repeats():
    seq1 = random_seq(1000)
    seq3 = seq1[500:520]
    seq2 = (
        seq1[:500] + seq3 + seq1[500 - len(seq3) : 600] + seq3 + seq1[600 + len(seq3) :]
    )
    assert complexity_score(seq2) > 10 + complexity_score(seq1)
