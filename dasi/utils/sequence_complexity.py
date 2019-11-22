"""Module for calculating the complexity of a DNASequence.

Complexity rules borrowed from IDT.

* less than 25% GC content
* higher than 75% GC content
* homopolymeric runs of 10 or more A/Ts
* homopolymeric runs of 6 or more G/Cs
* repeats greater than 14bp
* hairpins
* GC scew?
* windowed scew?
"""
import re
from typing import Dict
from typing import Union


class ComplexityConfig:
    HOMOPOLYMERIC_AT = 6
    HOMOPOLYMERIC_GC = 10
    REPEAT_MAX_LENGTH = 14
    PARTITION_OVERLAP = 25


C = ComplexityConfig


at_pattern = re.compile("[AT]{" + str(C.HOMOPOLYMERIC_AT) + ",}")
gc_pattern = re.compile("[GC]{" + str(C.HOMOPOLYMERIC_GC) + ",}")


def get_AT_complexity(seq):
    cmax = 0
    for match in at_pattern.findall(seq, re.IGNORECASE):
        if len(match) >= C.HOMOPOLYMERIC_AT:
            c = len(match) - C.HOMOPOLYMERIC_AT + 1
            if c > cmax:
                cmax = c
    return cmax


def get_GC_complexity(seq):
    cmax = 0
    for match in gc_pattern.findall(seq, re.IGNORECASE):
        if len(match) >= C.HOMOPOLYMERIC_GC:
            c = len(match) - C.HOMOPOLYMERIC_GC + 1
            if c > cmax:
                cmax = c
    return cmax


def get_GC_content(seq):
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def get_GC_content_complexity(seq):
    gc = get_GC_content(seq)
    return abs(gc * 100.0 - 50) * 17 / 25.0


def iter_kmers(seq, length):
    for i in range(0, len(seq) - length):
        yield (i, seq[i : i + length])


def count_kmers(seq, length):
    kmers = {}
    for kmer in iter_kmers(seq, length):
        kmers.setdefault(kmer[1], list())
        kmers[kmer[1]].append(kmer[0])
    return kmers


def get_repeat_complexity(seq, length):
    repeats = {}
    for kmer, positions in count_kmers(seq, length).items():
        if len(positions) > 1:
            kmer1 = kmer2 = kmer
            length = len(kmer)
            while kmer1 == kmer2:
                length += 1
                kmer1 = seq[positions[0] : positions[0] + length]
                kmer2 = seq[positions[1] : positions[1] + length]
            passes = True
            if len(kmer1) >= length:
                for k in repeats:
                    if kmer1 in k:
                        passes = False
                        break
                if passes:
                    repeats[kmer1] = (len(kmer1) - length) * 2
    return repeats


def complexity(seq: str) -> Dict[str, Union[float, int]]:
    complexity_info = {}
    complexity_info["Homopolymeric AT"] = get_AT_complexity(seq)
    complexity_info["Homopolymeric GC"] = get_GC_complexity(seq)
    complexity_info["High GC Content"] = get_GC_content_complexity(seq)
    complexity_info["Repeats"] = get_repeat_complexity(seq, C.REPEAT_MAX_LENGTH)
    return complexity_info


def complexity_score(seq: str) -> float:
    data = complexity(seq)
    total = 0
    for k, v in data.items():
        if isinstance(v, int) or isinstance(v, float):
            total += v
        else:
            for _v in v.values():
                total += _v
    return total
