import pytest

from dasi.utils.biopython import random_sequence
from dasi.utils.sequence import DNAStats
from dasi.utils.sequence.sequence_partitioner import find_best_partitions
from dasi.utils.sequence.sequence_partitioner import find_by_partitions_for_sequence
from dasi.utils.sequence.sequence_partitioner import find_fast_opt_partition
from dasi.utils.sequence.sequence_partitioner import find_opt_partition


@pytest.mark.parametrize("delta", [None, 10])
@pytest.mark.parametrize("step_size", [10])
def test_find_opt_partition(delta, step_size):
    repeat = random_sequence(30)

    seq = (
        random_sequence(1000)
        + repeat
        + random_sequence(20)
        + repeat
        + random_sequence(1000)
    )
    stats = DNAStats(seq, 20, 20, 20)

    p, cmin = find_opt_partition(stats, 10, step_size=step_size, delta=delta)
    assert p > 1000 and p < 1000 + 30 + 20 + 30


@pytest.mark.parametrize("delta", [None, 10])
@pytest.mark.parametrize("step_size", [10])
def test_find_fast_opt_partition(delta, step_size):
    repeat = random_sequence(30)

    seq = (
        random_sequence(1000)
        + repeat
        + random_sequence(20)
        + repeat
        + random_sequence(1000)
    )
    stats = DNAStats(seq, 20, 20, 20)

    p, _ = find_fast_opt_partition(
        stats, 10, threshold=10, step_size=step_size, delta=delta
    )
    assert p > 1000 and p < 1000 + 30 + 20 + 30


@pytest.mark.parametrize("delta", [None, 10])
@pytest.mark.parametrize("step_size", [10])
def test_find_best_partitions(delta, step_size):
    repeat = random_sequence(30)

    seq = (
        random_sequence(1000)
        + repeat
        + random_sequence(20)
        + repeat
        + random_sequence(1000)
    )
    stats = DNAStats(seq, 20, 20, 20)

    p = find_best_partitions(stats, threshold=10, step_size=step_size, delta=delta)
    for _p in p:
        assert _p > 1000 and _p < 1000 + 30 + 20 + 30


@pytest.mark.parametrize("delta", [None, 10])
@pytest.mark.parametrize("step_size", [10])
@pytest.mark.parametrize("cyclic", [True, False])
def test_find_partitions_for_sequence(delta, step_size, cyclic):
    repeat = random_sequence(30)

    seq = (
        random_sequence(1000)
        + repeat
        + random_sequence(20)
        + repeat
        + random_sequence(1000)
    )
    stats = DNAStats(seq, 20, 20, 20)

    p = find_by_partitions_for_sequence(
        stats, cyclic, threshold=10, step_size=step_size, delta=delta
    )
    passes = False
    for _p in p:
        if _p > 1000 and _p < 1000 + 30 + 20 + 30:
            passes = True
    assert passes


def test_():
    backbone = random_sequence(3000)
    repeat = random_sequence(30)
    complex_sequence = repeat + random_sequence(200) + repeat + random_sequence(1000)
    goal = backbone[1000:] + complex_sequence + backbone[:1000]

    stats = DNAStats(goal, 14, 20, 20)
    p = find_by_partitions_for_sequence(
        stats, cyclic=True, threshold=10, step_size=10, delta=30
    )
    print(p)
