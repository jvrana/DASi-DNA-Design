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

import numpy as np


class ComplexityConfig:
    HOMOPOLYMERIC_AT = 10
    HOMOPOLYMERIC_GC = 6
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
    return cmax / 10.0


def get_GC_complexity(seq):
    cmax = 0
    for match in gc_pattern.findall(seq, re.IGNORECASE):
        if len(match) >= C.HOMOPOLYMERIC_GC:
            c = len(match) - C.HOMOPOLYMERIC_GC + 1
            if c > cmax:
                cmax = c
    return cmax / 10.0


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


def repeat_complexity_func(len_repeat, max_repeat_length):
    return (len_repeat - max_repeat_length) / 2.0 + 10


def get_repeat_complexity(seq, length):
    repeats = {}
    for kmer, positions in count_kmers(seq, length).items():
        if len(positions) > 1:
            kmer1 = kmer2 = kmer
            kmer_len = len(kmer)
            while kmer1 == kmer2:
                kmer_len += 1
                kmer1 = seq[positions[0] : positions[0] + kmer_len]
                kmer2 = seq[positions[1] : positions[1] + kmer_len]
            passes = True
            if len(kmer1) >= length:
                for k in repeats:
                    if kmer1 in k:
                        passes = False
                        break
                if passes:
                    repeats[kmer1] = repeat_complexity_func(len(kmer1), length)
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


class DNAStats:

    BASES = "AGCT"

    def __init__(self, seq, repeat_window, stats_window, hairpin_window):
        self.seq = seq
        arr = np.array(list(seq))
        self.bases = self.BASES

        self.seq_onehot = self.one_hot(arr, self.bases)
        self.repeat_window = repeat_window
        self.stats_window = stats_window
        self.repeat_signatures = self.get_repeat_signatures(repeat_window)

        rolling_stats = self.get_base_stats(stats_window)
        gc_content = (
            rolling_stats[self.BASES.index("G"), :]
            + rolling_stats[self.BASES.index("C"), :]
        )
        at_content = (
            rolling_stats[self.BASES.index("A"), :]
            + rolling_stats[self.BASES.index("T"), :]
        )
        self.rolling_stats = np.vstack((rolling_stats, gc_content, at_content))

        self.fwd_signatures, self.rev_signatures = self.get_hairpin_signatures(
            hairpin_window
        )

    @staticmethod
    def one_hot(sequence, categories):
        arrs = []
        for i, c in enumerate(categories):
            a = sequence == c
            a = a.astype(int)
            arrs.append(a)
        return np.vstack(arrs)

    @staticmethod
    def rolling_window(a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

    @staticmethod
    def rolling_average(x, n):
        mv = np.convolve(x, np.ones(n) / n, mode="valid")
        return np.concatenate(([np.NaN for k in range(n - 1)], mv))

    def get_base_stats(self, window):
        return np.mean(self.rolling_window(self.seq_onehot, window), axis=2)

    def get_repeat_signatures(self, window):
        base_signatures = np.random.uniform(0.0, 1.0, size=4).reshape(-1, 4)
        seq_signatures = np.sum(np.multiply(self.seq_onehot, base_signatures.T), axis=0)
        x = np.random.uniform(0.0, 100.0, size=window)
        mv = np.convolve(seq_signatures, x, mode="valid")
        signatures = np.concatenate(([np.NaN for _ in range(window - 1)], mv))
        return signatures

    def get_hairpin_signatures(self, window):
        base_signatures = np.random.uniform(0.0, 1.0, size=4).reshape(-1, 4)
        seq_signatures = np.sum(np.multiply(self.seq_onehot, base_signatures.T), axis=0)
        revcomp_signatures = np.sum(
            np.multiply(self.seq_onehot[:, ::-1], base_signatures[:, ::-1].T), axis=0
        )
        x = np.random.uniform(0.0, 100.0, size=window)
        mv1 = np.convolve(seq_signatures, x, mode="valid")
        mv2 = np.convolve(revcomp_signatures, x, mode="valid")
        return mv1, mv2[::-1]

    @staticmethod
    def count_signatures(a):
        return len(np.where(np.unique(a, return_counts=True)[1] > 1)[0])

    def slice_repeats(self, i, j):
        return self.count_signatures(self.repeat_signatures[i:j])

    def slice_hairpins(self, i, j):
        f = self.fwd_signatures[i:j]
        r = self.rev_signatures[i:j]
        d = np.concatenate((f, r))
        num_hairpins = self.count_signatures(d) - 2 * self.count_signatures(f)
        return num_hairpins

    def window_cost(self, i, j):
        d = self.rolling_stats[:4, i:j]

        # keep each base below 0.8
        a = np.sum(d > 0.75, axis=1)

        # keep gc content below extremes 0.7
        b = np.sum(self.rolling_stats[4, i:j] > 0.85)
        return a, b

    def __call__(self, i=None, j=None):
        return {
            "n_repeats": self.slice_repeats(i, j),
            "n_hairpins": self.slice_hairpins(i, j),
            "window_cost": self.window_cost(i, j),
        }

    def cost(self, i=None, j=None):
        d = self(i, j)
        w1 = d["window_cost"][0]
        w2 = d["window_cost"][1]
        total = d["n_repeats"] + d["n_hairpins"] + np.sum(w1) + np.sum(w2)
        return total
