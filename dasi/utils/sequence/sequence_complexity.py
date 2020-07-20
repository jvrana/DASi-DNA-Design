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
from functools import lru_cache

import numpy as np


# TODO:' compute cost for extremes of GC content
# TODO: this could be a separate library
class DNAStats:
    """
    Class for computing statistics for DNA Sequences.

    ::

        {
            "n_repeats",
            "n_hairpins",
            "window_cost",
            "gc_cost"
        }
    """

    BASES = "AGCT"
    np.random.seed(1)
    BASE_SIGNATURES = np.random.uniform(0.0, 1.0, size=len(BASES)).reshape(-1, 4)

    DEFAULT_MODE = ("window", "hairpin", "repeat")
    ONLY_WINDOW = ("window",)
    ONLY_HAIRPIN = ("hairpin",)
    ONLY_REPEAT = "repeat"

    def __init__(
        self,
        seq,
        repeat_window,
        stats_window,
        hairpin_window,
        base_percentage_threshold: float = 0.75,
        gc_content_threshold: float = 0.85,
        at_content_threshold: float = 0.90,
        conv_seed_repeat: np.array = None,
        conv_seed_hairpin: np.array = None,
        mode: str = None,
    ):
        """

        :param seq:
        :param repeat_window: the length of the k-mers to search for repeats.
        :param stats_window: the width of the stats window
        :param hairpin_window: length of hte k-mers to search for hairpins.
        :param base_percentage_threshold: if the per base percentage threshold
            to count in slinding window calculatiions. For example,
            for a sliding window of 20, if the threshold is 0.75, if 15 bases
            or of a single base pair, this will get counted in the `window_cost`
            method.
        :param gc_content_threshold: if the gc percentage threshold
            to count in sliding window calculatiions. For example,
            for a sliding window of 20, if the threshold is 0.75, if 15 bases
            of any window of 20 were gc, this will get counted in the `window_cost`
            method.
        :param conv_seed_repeat: (optional) the array to convolve the one-hot
            encoded DNA sequence with. If not provided, will be generated by
            `np.random.uniform(0.0, 100.0, size=repeat_window)`
        :param conv_seed_hairpin: (optional) the array to convolve the one-hot
            encoded DNA sequence with to find hairpin signatures.
            If not provided, will be generated by
            `np.random.uniform(0.0, 100.0, size=hairpin_window)`
        :param mode: (optional) Set to DNAStats.ONLY_HAIRPIN to only
            calculate hairpin signatures, DNAStats.ONLY_REPEATS to only calculate
            repeat signatures. DNAStats.ONLY_WINDOW to perform only
            sliding window calculations.
        """
        assert isinstance(seq, str)
        seq = seq.upper()
        self.seq = seq
        arr = np.array(list(self.seq))
        self.bases = self.BASES.upper()

        self.seq_onehot = self.one_hot(arr, self.bases)
        self.repeat_window = repeat_window
        self.stats_window = stats_window
        self.hairpin_window = hairpin_window
        self.gc_content_threshold = gc_content_threshold
        self.at_content_threshold = at_content_threshold
        self.base_percentage_threshold = base_percentage_threshold
        if not mode:
            mode = self.DEFAULT_MODE
        self.mode = mode
        if conv_seed_repeat:
            self.conv_seed_repeat = conv_seed_repeat
        else:
            self.conv_seed_repeat = np.random.uniform(0.0, 100.0, size=repeat_window)
        if conv_seed_hairpin:
            self.conv_seed_hairpin = conv_seed_hairpin
        else:
            self.conv_seed_hairpin = np.random.uniform(0.0, 100.0, size=hairpin_window)

        if "window" in self.mode:
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
        else:
            self.rolling_stats = None

        ###
        # get repeat signatures
        ###
        seq_signatures = np.sum(
            np.multiply(self.seq_onehot, self.BASE_SIGNATURES.T), axis=0
        )
        if "repeat" in self.mode:
            mv = np.convolve(seq_signatures, self.conv_seed_repeat, mode="valid")
            self.repeat_signatures = np.concatenate(
                ([np.NaN for _ in range(repeat_window - 1)], mv)
            )
        else:
            self.repeat_signatures = None

        ###
        # get hairpin signatures
        ###
        if "hairpin" in self.mode:
            revcomp_signatures = np.sum(
                np.multiply(self.seq_onehot[:, ::-1], self.BASE_SIGNATURES[:, ::-1].T),
                axis=0,
            )
            self.fwd_signatures = np.convolve(
                seq_signatures, self.conv_seed_hairpin, mode="valid"
            )
            self.rev_signatures = np.convolve(
                revcomp_signatures, self.conv_seed_hairpin, mode="valid"
            )[::-1]
        else:
            self.fwd_signatures = None
            self.rev_signatures = None

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

    def window_cost(
        self,
        i,
        j,
        base_threshold_perc: float = None,
        gc_content_threshold: float = None,
        at_content_threshold: float = None,
    ):

        if base_threshold_perc is None:
            base_threshold_perc = self.base_percentage_threshold
        if gc_content_threshold is None:
            gc_content_threshold = self.gc_content_threshold

        if at_content_threshold is None:
            at_content_threshold = self.at_content_threshold
        d = self.rolling_stats[:4, i:j]

        # keep each base below 0.8
        a = np.sum(d > base_threshold_perc, axis=1)

        # keep gc content below extremes 0.7
        b = np.sum(self.rolling_stats[4, i:j] > gc_content_threshold)

        # keep gc content below extremes 0.7
        c = np.sum(self.rolling_stats[5, i:j] > at_content_threshold)
        return a, b, c

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

    def count_repeats_from_slice(self, i, j):
        """Counts the number of repeats found in slices [0, i) and [j, None) of
        all kmers in the sequence found within the [i, j) slice.

        .. note::
            The hairpin_window_length of the DNAStats instance is used as the kmer
            length.

        :param i: start of the slice (inclusive)
        :param j: end of the slice (exclusive)
        :return: tuple of repeat counts found in left and right slices respectively.
        """
        if i < 0:
            raise IndexError("i={} cannot be < 0")
        if j > len(self.seq):
            raise IndexError("j={} cannot be >= length {}".format(j, len(self.seq)))
        kmer_length = self.hairpin_window

        # we ust signatures (src) of kmers inside i:j...
        src_indices = np.arange(i, j - kmer_length + 1)

        # too look for signatures (dest) in the left and right flanks
        left_flank_indices = np.arange(0, i)
        right_flank_indices = np.arange(j, len(self.seq) - kmer_length)
        dest_indices = np.concatenate([left_flank_indices, right_flank_indices])

        src = self.fwd_signatures[src_indices]

        # TODO: unsure why the rev indices need to be shifted...
        src_rc = self.rev_signatures[src_indices]
        dest = self.fwd_signatures[dest_indices]

        def _mispriming_counts(src_sig, src_sig_rc, dest_sig):
            # 5' end, top repeats to self
            counts_trim_top = self.count_signatures(src_sig)

            # 5' end, calculate priming to top strand
            x, counts = np.unique(src_sig, return_counts=True)
            top_counts = counts[np.where(np.isin(x, dest_sig))].sum()

            # 5' end, calculate priming to bottom strand
            x, counts = np.unique(src_sig_rc, return_counts=True)
            bottom_counts = (
                counts[np.where(np.isin(x, dest_sig))].sum() + (counts - 1).sum()
            )

            return top_counts + bottom_counts + counts_trim_top

        return _mispriming_counts(src, src_rc, dest)

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


@lru_cache(512)
def cached_stats(seq: str, hairpin_window: int):
    return DNAStats(seq, 1, 1, hairpin_window, mode=DNAStats.ONLY_HAIRPIN)


@lru_cache(512)
def count_misprimings_in_amplicon(
    seq,
    i,
    j,
    min_primer_anneal: int = 12,
    max_primer_anneal: int = 30,
    cyclic: bool = False,
    max_length: int = 3000,
):
    """Counts the estimated number of misprimings in an amplicon.

    :param seq: The template sequence
    :param i: amplicon start (inclusive)
    :param j: amplicon end (exclusive)
    :param min_primer_anneal: The minimum number of bases to look for misprimings
        (exact matches only).
    :param max_primer_anneal: The maximum size of the annealing sequence of the primers.
        Internally, this takes a window slice from the ends of the amplicon to look
        for repeats in the template.
    :param cyclic: If True, assumes the sequence is cyclic.
    :param max_length: The maximum length to consider a mispriming. Generally,
        this should be about the bp distinguishable in your fragment purification
        (e.g. gel)
    :return:
    """
    if j == i:
        return 0
    if i < 0:
        raise IndexError("i cannot be < 0")
    if j < 0:
        raise IndexError("j cannot be < 0")
    elif j < i and not cyclic:
        raise IndexError(
            "Invalid indices provided. "
            "Indices [{i},j) indicates a cyclic sequence, but "
            "cyclic was set to False."
        )
    elif cyclic:
        seq = seq[i:] + seq[:i]
        if j < i:
            length = len(seq) - i + j
            i = 0
            j = length
        else:
            j = j - i
            i = 0

        # TODO: trim extra long sequences
        # trim extra long sequences
        # if len(seq) - j > 2 * max_length:
        #     seq = seq[:j] + seq[j:][max_length] + seq[-max_length:]
    # else:
    #     # trim extra long sequences
    #     if i > max_length:
    #         delta = max_length - i
    #     else:
    #         delta = 0
    #     seq = seq[delta:i] + seq[i:j] + seq[j:][:max_length]
    #     if i > max_length:
    #         delta = max_length - i
    #         i -= delta
    #         j -= delta
    stats = cached_stats(seq, min_primer_anneal)
    n1 = stats.count_repeats_from_slice(i, min(i + max_primer_anneal, len(seq) - 1))
    k = j - max_primer_anneal
    if k < 0:
        n2 = 0
    else:
        n2 = stats.count_repeats_from_slice(max(j - max_primer_anneal, 0), j)
    return n1 + n2
