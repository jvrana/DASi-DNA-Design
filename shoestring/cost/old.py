import itertools
import json
import os
import time
from tqdm import tqdm
import numpy as np
from shoestring.log import logger


HERE = os.path.dirname(os.path.abspath(__file__))


def digitize_dictionary(mydict, min, max):
    bins = list(mydict.keys())
    x = np.arange(min, max + 1)
    key_indices = np.digitize(x, bins)
    values = []
    for index in key_indices:
        if index == len(bins):
            index = index - 1
        values.append(mydict[bins[index]])
    return dict(zip(x, values))


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        print(f"Starting {method.__name__}")
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print(f"{method.__name__} took {(te-ts)*1000} ms\n")
        return result

    return timed


class GibsonAssemblyCost(object):
    # -----------------------------------------------
    # Parameters
    # -----------------------------------------------

    PRIMER_EXTENSION = (0, 40)
    PRIMER_COST = 0.60  # / bp
    ULTRAMER_EXTENSION = (40, 180)
    ULTRAMER_COST = 0.80  # / bp
    INF = 1_000_000.0

    GENE_SYNTHESIS_COST_RANGES = [
        {"lower": 0, "upper": 0, "cost": 0.0},
        {"lower": 100, "upper": 500, "cost": 89.0},
        {"lower": 501, "upper": 750, "cost": 129.0},
        {"lower": 751, "upper": 1000, "cost": 149.0},
        {"lower": 1001, "upper": 1250, "cost": 209.0},
        {"lower": 1251, "upper": 1500, "cost": 249.0},
        {"lower": 1501, "upper": 1750, "cost": 289.0},
        {"lower": 1751, "upper": 2000, "cost": 329.0},
        {"lower": 2001, "upper": INF, "cost": INF},
    ]

    SYNTHESIS_GAP_SPAN_RANGE = (-300, 1000, 10)
    SYNTHESIS_SIZE_OPTIONS = (100, 2000, 10)

    JUNCTION_EFFICIENCY = [
        {"lower": -INF, "upper": 10, "cost": 0.0},
        {"lower": 11, "upper": 19, "cost": 0.1},
        {"lower": 20, "upper": 30, "cost": 0.75},
        {"lower": 31, "upper": 100, "cost": 0.9},
        {"lower": 101, "upper": INF, "cost": 0.0},
    ]

    JUNCTION_GAP_SPAN_RANGE = (-300, 1000, 5)

    def __init__(self):
        self._gap_cost_dict_both_extendable = None
        self._gap_cost_dict_one_extendable = None
        self._gap_cost_dict_none_extendable = None
        self._gap_cost_dict_synthesis_none_extendable = None
        self._gap_cost_dict_synthesis_one_extendable = None
        self._gap_cost_dict_synthesis_extendable = None

    # -----------------------------------------------
    # Properties
    # -----------------------------------------------
    @property
    def synthesis_gap_bins(self):
        return np.arange(*self.SYNTHESIS_GAP_SPAN_RANGE)

    @property
    def junction_gap_bins(self):
        return np.arange(*self.JUNCTION_GAP_SPAN_RANGE)

    # -----------------------------------------------
    # Primer, Ultramer, Synthesis options
    # -----------------------------------------------

    def primer_extension_options(self):
        primer_options = []
        for extension in range(*self.PRIMER_EXTENSION):
            primer_options.append(
                {
                    "extension": extension,
                    "cost": self.PRIMER_COST * (20.0 + extension),
                    "name": "primer",
                }
            )
        return primer_options
        # primer_options.append((0, 0, "no primer"))

    def ultramer_extension_options(self):
        ultramer_options = []
        for extension in range(*self.ULTRAMER_EXTENSION):
            ultramer_options.append(
                {
                    "extension": extension,
                    "cost": self.ULTRAMER_COST * (20.0 + extension),
                    "name": "ultra-mer",
                }
            )
        return ultramer_options

    def synthesis_options(self):
        synthesis_options = []
        for gene_size in range(*self.SYNTHESIS_SIZE_OPTIONS):
            synthesis_options.append(
                {
                    "size": gene_size,
                    "cost": self._synthesis_cost(gene_size),
                    "name": "synthesis",
                }
            )
        return synthesis_options

    # -----------------------------------------------
    # Methods
    # -----------------------------------------------

    @staticmethod
    def find_in_range(
        x, ranges, lower_label="lower", upper_label="upper", label="cost"
    ):
        for range in ranges:
            lower = range[lower_label]
            upper = range[upper_label]
            if lower <= x <= upper:
                return range[label]

    def _synthesis_cost(self, size):
        return self.find_in_range(size, self.GENE_SYNTHESIS_COST_RANGES)

    def _junction_efficiency(self, overlap):
        return self.find_in_range(overlap, self.JUNCTION_EFFICIENCY)

    def _optimize_primer(self, primer_choices):
        """
        Create a dictionary of gap_span to best_extension case.

        :param gap_span_start:
        :type gap_span_start:
        :param gap_span_end:
        :type gap_span_end:
        :param gap_span_step:
        :type gap_span_step:
        :param primer_choices:
        :type primer_choices:
        :return:
        :rtype:
        """
        results = {}

        # iterate through gap span
        for gap_span in tqdm(
            range(*self.JUNCTION_GAP_SPAN_RANGE), desc="optimizing primers"
        ):
            results[gap_span] = None
            for choice in primer_choices:
                if not isinstance(choice, tuple):
                    choice = (choice,)
                extensions = [x["extension"] for x in choice]
                costs = [x["cost"] for x in choice]
                names = [x["name"] for x in choice]

                total_extension = sum(extensions)
                total_cost = sum(costs)
                junction_efficiency = self._junction_efficiency(
                    -gap_span + total_extension
                )

                # score
                score = self.INF
                if junction_efficiency > 0:
                    score = total_cost / junction_efficiency

                _x = results[gap_span]
                if _x is None or score < _x["score"]:
                    results[gap_span] = {
                        "score": score,
                        "cost": total_cost,
                        "efficiency": junction_efficiency,
                        "condition": f"{names}",
                        "extension": f"{extensions}",
                    }
        return results

    def _optimize_synthesis_options(
        self, synthesis_options, left_ext_dict, right_ext_dict
    ):
        """

        :param gap_span_start:
        :type gap_span_start:
        :param gap_span_end:
        :type gap_span_end:
        :param gap_span_step:
        :type gap_span_step:
        :param synthesis_options:
        :type synthesis_options:
        :param extension_dict:
        :type extension_dict:
        :return:
        :rtype:
        """

        # for gap in gap_span_options
        # for gene in gene_options
        # for offset in offset_options
        # best primer option for left_offset
        # best primer option for right offset

        results = {}

        # digitize dictionary so that it includes all possible values
        largest_overlap = (
            self.SYNTHESIS_SIZE_OPTIONS[1] - self.SYNTHESIS_GAP_SPAN_RANGE[0]
        )
        smallest_overlap = (
            self.SYNTHESIS_SIZE_OPTIONS[0] - self.SYNTHESIS_GAP_SPAN_RANGE[1]
        )
        left_ext_dict = digitize_dictionary(
            left_ext_dict, -largest_overlap, -smallest_overlap
        )
        right_ext_dict = digitize_dictionary(
            right_ext_dict, -largest_overlap, -smallest_overlap
        )

        # optimize for gap span
        for gap_span in tqdm(
            range(*self.SYNTHESIS_GAP_SPAN_RANGE), desc="optimizing synthesis"
        ):
            results[gap_span] = None

            # optimize for gene size
            for gene in synthesis_options:
                max_overlap = gene["size"] - gap_span

                extension_params = {
                    "score": self.INF,
                    "cost": self.INF,
                    "eff": 0.0,
                    "left": None,
                    "right": None,
                    "gene_offset": None,
                }

                # optimize for offset
                for offset in range(0, max_overlap):
                    gap_left = -max_overlap + offset
                    gap_right = -max_overlap - gap_left
                    ext_left = left_ext_dict[gap_left]
                    ext_right = right_ext_dict[gap_right]
                    ext_cost = ext_left["cost"] + ext_right["cost"]
                    ext_eff = ext_left["efficiency"] * ext_right["efficiency"]
                    ext_score = self.INF
                    if ext_eff > 0:
                        ext_score = ext_cost / ext_eff

                    if ext_score < extension_params["score"]:
                        extension_params["score"] = ext_score
                        extension_params["cost"] = ext_cost
                        extension_params["eff"] = ext_eff
                        extension_params["left"] = ext_left
                        extension_params["right"] = ext_right
                        extension_params["gene_offset"] = offset

                cost = gene["cost"] + extension_params["cost"]
                eff = extension_params["eff"]
                score = self.INF
                if eff > 0:
                    score = cost / eff

                if results[gap_span] is None or score < results[gap_span]["score"]:
                    results[gap_span] = {
                        "score": score,
                        "cost": cost,
                        "gene_cost": gene["cost"],
                        "extension": extension_params,
                        "gene_size": gene["size"],
                        "eff": eff,
                    }
        return results

    def gap_cost_dict_both_extendable(self):
        if self._gap_cost_dict_both_extendable is None:
            both_extendable_options = list(
                itertools.combinations(
                    self.primer_extension_options() + self.ultramer_extension_options(),
                    2,
                )
            )
            self._gap_cost_dict_both_extendable = self._optimize_primer(
                both_extendable_options
            )
        return self._gap_cost_dict_both_extendable

    def gap_cost_dict_one_extendable(self):
        if self._gap_cost_dict_one_extendable is None:
            self._gap_cost_dict_one_extendable = self._optimize_primer(
                self.primer_extension_options() + self.ultramer_extension_options()
            )
        return self._gap_cost_dict_one_extendable

    def gap_cost_dict_none_extendable(self):
        if self._gap_cost_dict_none_extendable is None:
            self._gap_cost_dict_none_extendable = self._optimize_primer(
                [{"extension": 0, "cost": 1.0, "name": "no primer"}]
            )
        return self._gap_cost_dict_none_extendable

    def _gap_cost_dict_synthesis(self, left_ext=False, right_ext=False):
        left = right = self.gap_cost_dict_none_extendable()
        if left_ext:
            left = self.gap_cost_dict_one_extendable()
        if right_ext:
            right = self.gap_cost_dict_one_extendable()
        return self._optimize_synthesis_options(self.synthesis_options(), left, right)

    def gap_cost_dict_synthesis_none_extendable(self):
        if self._gap_cost_dict_synthesis_none_extendable is None:
            self._gap_cost_dict_synthesis_none_extendable = self._gap_cost_dict_synthesis(
                left_ext=False, right_ext=False
            )
        return self._gap_cost_dict_synthesis_none_extendable

    def gap_cost_dict_synthesis_one_extendable(self):
        if self._gap_cost_dict_synthesis_one_extendable is None:
            self._gap_cost_dict_synthesis_one_extendable = self._gap_cost_dict_synthesis(
                left_ext=True, right_ext=False
            )
        return self._gap_cost_dict_synthesis_one_extendable

    def gap_cost_dict_synthesis_extendable(self):
        if self._gap_cost_dict_synthesis_extendable is None:
            self._gap_cost_dict_synthesis_extendable = self._gap_cost_dict_synthesis(
                left_ext=True, right_ext=True
            )
        return self._gap_cost_dict_synthesis_extendable

    @classmethod
    def _path_name(cls):
        name = cls.__name__
        path = os.path.join(HERE, name + ".json")
        return path

    def save(self):
        with open(self._path_name(), "w") as handle:
            json.dump(self.__dict__, handle)

    @classmethod
    def load(cls):
        if os.path.isfile(cls._path_name()):
            with open(cls._path_name(), "rb") as handle:
                data = json.load(handle)
            x = cls()
            x.__dict__.update(data)
            return x

    def compute(self):

        logger.debug("Computing gap cost none extendable")
        self.gap_cost_dict_none_extendable()

        logger.debug("Computing gap cost one extendable")
        self.gap_cost_dict_one_extendable()

        logger.debug("Computing gap cost both extendable")
        self.gap_cost_dict_both_extendable()

        logger.debug("Computing synthesis cost none extendable")
        self.gap_cost_dict_synthesis_none_extendable()

        logger.debug("Computing synthesis cost one extendable")
        self.gap_cost_dict_synthesis_one_extendable()

        logger.debug("Computing synthesis cost both extendable")
        self.gap_cost_dict_synthesis_extendable()

    def conditions(self):
        dicts = [
            [0, False, self.gap_cost_dict_none_extendable()],
            [1, False, self.gap_cost_dict_one_extendable()],
            [2, False, self.gap_cost_dict_both_extendable()],
            [0, True, self.gap_cost_dict_synthesis_none_extendable()],
            [1, True, self.gap_cost_dict_synthesis_one_extendable()],
            [2, True, self.gap_cost_dict_synthesis_extendable()],
        ]

        for d in dicts:
            d[-1] = digitize_dictionary(d[-1], -500, 5000)
        return dicts

    def filter_by_condition(self, n, syn):
        dicts = self.conditions()
        if not syn:
            dicts = list(filter(lambda x: not x[1], dicts))
        return list(filter(lambda x: x[0] == n, dicts))

    def gap_cost_dict(self, min, max, e=0, syn=True):
        dicts = self.filter_by_condition(e, syn)
        data = {}
        for span in tqdm(range(min, max + 1)):
            values = np.array([d[-1][span]["score"] for d in dicts])
            data[span] = np.min(values)
        return data

    def gap_cost(self, span, e=0, syn=True):
        dicts = self.filter_by_condition(e, syn)

        values = [d[-1][span]["score"] for d in dicts]
        index = values.index(min(values))
        return dicts[index][-1][span]["score"]

        # condition 1: no synthesis, both extendable
        # condition 2: no synthesis, one extendable
        # condition 3: no synthesis, none extendable
        # condition 4: synthesis, both extendable
        # condition 5: synthesis, one extendable
        # condition 6: synthesis, none extendable


# gac = GibsonAssemblyCost()
# print("loading")
# gac.load()
# print('computing')
# gac.compute()
# gac.save()

# d = gac.gap_cost_dict(-200, 1000)
# pass
#
#
# # gac.compute()
# # gac.save()
# #
# r1 = gac.gap_cost_dict_synthesis_none_extendable()
# for k, v in r1.items():
#     print(f"{k} {v}")
#
# r2 = gac.gap_cost_dict_synthesis_one_extendable()
# for k, v in r2.items():
#     print(f"{k} {v}")
#
# r3 = gac.gap_cost_dict_synthesis_extendable()
# for k, v in r3.items():
#     print(f"{k} {v}")
#
# r4 = gac.gap_cost_dict_both_extendable()
# for k, v in r4.items():
#     print(f"{k} {v}")
#
# r5 = gac.gap_cost_dict_one_extendable()
# for k, v in r5.items():
#     print(f"{k} {v}")
#
# r6 = gac.gap_cost_dict_none_extendable()
# for k, v in r6.items():
#     print(f"{k} {v}")
#
# import matplotlib.pyplot as plt
#
# X = list(range(-400, 1000, 5))
# Y = []
# for i in X:
#     Y.append(gac.gap_cost(i)['score'])
#
#
# plt.figure()
# plt.scatter(X, Y, c='brown')
# plt.scatter(list(r1.keys()), [v['score'] for _, v in r1.items()], c='r')
# plt.scatter(list(r2.keys()), [v['score'] for _, v in r2.items()], c='b')
# plt.scatter(list(r3.keys()), [v['score'] for _, v in r3.items()], c='g')
# plt.scatter(list(r4.keys()), [v['score'] for _, v in r4.items()], c='c')
# plt.scatter(list(r5.keys()), [v['score'] for _, v in r5.items()], c='purple')
# plt.scatter(list(r6.keys()), [v['score'] for _, v in r6.items()], c='yellow')
# plt.ylim(0, 500)
# plt.show()
