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


# TODO: compute cost for extremes of GC content
class DNAStats:

    BASES = "AGCT"
    np.random.seed(1)
    BASE_SIGNATURES = np.random.uniform(0.0, 1.0, size=len(BASES)).reshape(-1, 4)

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
        self.cached_partitions = []

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
        seq_signatures = np.sum(
            np.multiply(self.seq_onehot, self.BASE_SIGNATURES.T), axis=0
        )
        x = np.random.uniform(0.0, 100.0, size=window)
        mv = np.convolve(seq_signatures, x, mode="valid")
        signatures = np.concatenate(([np.NaN for _ in range(window - 1)], mv))
        return signatures

    def get_hairpin_signatures(self, window):
        seq_signatures = np.sum(
            np.multiply(self.seq_onehot, self.BASE_SIGNATURES.T), axis=0
        )
        revcomp_signatures = np.sum(
            np.multiply(self.seq_onehot[:, ::-1], self.BASE_SIGNATURES[:, ::-1].T),
            axis=0,
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

    # def get_GC_content_complexity(seq):
    #     gc = get_GC_content(seq)
    #     return abs(gc * 100.0 - 50) * 17 / 25.0

    @staticmethod
    def _optimize_partition_helper(
        signatures: np.ndarray, step: int, i: int = None, j: int = None
    ):
        """Optimize partition by minimizing the number of signatures in the
        given array.

        :param signatures: array of signatures
        :param step: step size
        :param i:
        :param j:
        :return:
        """
        d = []

        if i is None:
            i = 0
        if j is None:
            j = signatures.shape[1]

        for x in range(i, j, step):
            m1 = np.empty(signatures.shape[1])
            m2 = m1.copy()
            m1.fill(np.nan)
            m2.fill(np.nan)

            m1[:x] = np.random.uniform(1, 10)
            m2[x:] = np.random.uniform(1, 10)

            d += [m1, m2]
        d = np.vstack(d)
        z = np.tile(d, signatures.shape[0]) * signatures.flatten()

        partition_index = np.repeat(
            np.arange(0, signatures.shape[1], step),
            signatures.shape[0] * signatures.shape[1] * 2,
        )

        a, b, c = np.unique(z, return_counts=True, return_index=True)
        i = b[np.where(c > 1)]
        a, c = np.unique(partition_index[i], return_counts=True)
        if len(c):
            arg = c.argmin()
            return a[arg], c[arg]

    def partition(
        self,
        step,
        overlap,
        i=None,
        j=None,
        border=100,
        stopping_threshold=10,
        window=None,
        use_cache=True,
    ):
        if j is None:
            j = len(self.seq)
        if i is None:
            i = 0

        if use_cache:
            for cached in self.cached_partitions:
                if i >= cached[0] and j <= cached[1]:
                    if cached[2]["cost"] < stopping_threshold:
                        window = (
                            cached[2]["index_1"][1] - overlap,
                            cached[2]["index_2"][0] + overlap,
                        )
                elif i <= cached[0] and j >= cached[1]:
                    if cached[2]["cost"] > stopping_threshold:
                        return []

        if not window:
            window = (i, j)

        to_search = []
        for index in range(window[0], window[1], step):
            i1 = index
            i2 = index - overlap
            if i2 <= border:
                continue
            if i1 > len(self.seq) - border:
                continue
            to_search.append((i1, i2))

        if not to_search:
            return []

        search_index = int(len(to_search) / 2)
        searched_costs = []
        searched = []
        while search_index not in searched:

            searched.append(search_index)

            try:
                i1, i2 = to_search[search_index]
            except IndexError:
                print(search_index)
                print(len(to_search))
                raise IndexError

            c1 = self.cost(None, i1)
            c2 = self.cost(i2, None)
            c = c1 + c2

            searched_costs.append((c, c1, c2, i1, i2, search_index))

            if c1 > c2:
                search_index = int(search_index / 2)
            elif c2 > c1:
                search_index = int((len(to_search) - search_index) / 2) + search_index

            if c < stopping_threshold:
                break

        searched_costs = sorted(searched_costs)
        c, c1, c2, i1, i2, search_index = searched_costs[0]

        partitions = self.partition_scan(
            1,
            overlap=overlap,
            i=i,
            j=j,
            window=(min(i1, i2) - overlap, max(i1, i2) + overlap),
        )
        if use_cache:
            self.cached_partitions.append((i, j, partitions[0]))
        return partitions

    def partition_scan(self, step, overlap, i=None, j=None, window=None):
        if j is None:
            j = len(self.seq)
        if i is None:
            i = 0

        if not window:
            window = (i, j)

        costs = []
        for x in range(window[0], window[1], step):
            i2, j2 = x, x - overlap

            c1 = self.cost(i, i2)
            c2 = self.cost(j2, j)
            costs.append(
                {
                    "cost": c1 + c2,
                    "cost_1": c1,
                    "cost_2": c2,
                    "index_1": (i, i2),
                    "index_2": (j2, j),
                }
            )
        costs = sorted(costs, key=lambda x: x["cost"])
        return costs

    def cost(self, i=None, j=None):
        d = self(i, j)
        w1 = d["window_cost"][0]
        w2 = d["window_cost"][1]
        total = (
            d["n_repeats"] + d["n_hairpins"] + np.sum(w1) + np.sum(w2) + d["gc_cost"]
        )
        return total

    def __call__(self, i=None, j=None):
        mn = np.mean(self.seq_onehot[:, i:j], axis=1)
        gi = self.bases.index("G")
        ci = self.bases.index("C")
        gc = mn[gi] + mn[ci]
        gccost = abs(gc * 100.0 - 50) * 17 / 25.0
        return {
            "n_repeats": self.slice_repeats(i, j),
            "n_hairpins": self.slice_hairpins(i, j),
            "window_cost": self.window_cost(i, j),
            "gc_cost": gccost,
        }
