"""sequence.

Methods for DNA/RNA sequence manipulation.
"""
from typing import List
from typing import Tuple

import numpy as np
from primer3plus.utils import anneal as p3panneal

BASES = "AGCT"


def random_sequence(size: int, letters=BASES) -> str:
    """Generate a random sequence string."""
    NP_BASES = np.array(list(letters))
    indices = np.random.randint(0, 4, size=size)
    bases = NP_BASES[indices]
    return "".join(bases)


def c(seq: str):
    """Return complement a dna string."""
    d = dict(zip("AGCTagct", "TCGAtcga"))
    return "".join([d[x] for x in seq])


def rc(seq: str):
    """Return a reverse complement a dna string."""
    return c(seq)[::-1]


def dna_like(seqstr: str, letters: str = "AGTCagctnNuU", min_length: int = 10):
    if seqstr is None:
        return False
    if len(seqstr) <= min_length:
        return False
    for s in seqstr:
        if s not in letters:
            return False
    return True


def anneal(
    template: str, primers: List[str], ignore_case: bool = True
) -> Tuple[List[dict], List[dict]]:
    assert isinstance(template, str)
    for p in primers:
        assert isinstance(p, str)
    return p3panneal(template, primers, ignore_case=ignore_case)
