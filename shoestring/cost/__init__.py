# TODO: convert everything to numpy math

from shoestring.utils import SpanDictionary, InfiniteDict
from tqdm import tqdm
from itertools import combinations
class Cost(object):
    INF = 100_000.0

    SCORING = {"materials": 1.0, "time": 20.0, "efficiency": 1.0}

    MIN_PRIMER_ANNEAL = 20

    PRIMER_COST_PER_BP = SpanDictionary(
        {
            (0, 40): {"base": 0, "bp": 0.6, "time": 1.5},
            (40, 80): {"base": 0, "bp": 0.8, "time": 2.0},
        }
    )

    PRIMER_COST = 0.60

    GENE_SYNTHESIS_COST_RANGES = SpanDictionary(
        {
            (0, 1): {"base": 0.0, "time": 0},
            (1, 100): {"base": INF, "time": INF},
            (100, 500): {"base": 89.0, "time": 3.0},
            (500, 750): {"base": 129.0, "time": 3.0},
            (750, 1000): {"base": 149.0, "time": 4.0},
            (1000, 1250): {"base": 209.0, "time": 7.0},
            (1250, 1500): {"base": 249.0, "time": 7.0},
            (1500, 1750): {"base": 289.0, "time": 7.0},
            (1750, 2000): {"base": 329.0, "time": 7.0},
            (2000, 2250): {"base": 399.0, "time": 7.0},
        }
    )

    JUNCTION_EFFICIENCY = SpanDictionary(
        {
            (-INF, 10): 0.0,
            (10, 20): 0.1,
            (20, 30): 0.75,
            (30, 100): 0.9,
            (100, INF): 0.0,
        }
    )

    JUNCTION_GAP_SPAN_RANGE = (-300, 1000, 5)

    def primer_extension_options(self):
        return dict(self.PRIMER_COST_PER_BP.items())

    @classmethod
    def score(cls, materials, time, efficiency):
        return {"materials": materials, "time": time, "efficiency": efficiency}

    def _primer_cost(self, v, primer_length):
        """Compute the cost of the primer from the parameter dictionaries above."""
        return v["base"] + v["bp"] * primer_length

    def score_primer_extension_cost(
        self, span, primer1=None, ext1=None, primer2=None, ext2=None
    ):
        """

        :param span:
        :param primer1: {'base', 'time', 'bp'}
        :param primer2: {'base', 'time', 'bp'}
        :param ext1:
        :param ext2:
        :return:
        """
        if primer1:
            primer1_cost = self._primer_cost(primer1, self.MIN_PRIMER_ANNEAL + ext1)
            t1 = primer1["time"]
        else:
            primer1_cost = 0
            ext1 = 0
            t1 = 0
        if primer2:
            primer2_cost = self._primer_cost(primer2, self.MIN_PRIMER_ANNEAL + ext2)
            t2 = primer2["time"]
        else:
            primer2_cost = 0
            ext2 = 0
            t2 = 0
        junction_efficiency = self.JUNCTION_EFFICIENCY[-span + ext1 + ext2][0]

        m = max([5.0, primer1_cost + primer2_cost])
        t = max([0.1, t1, t2])
        score_dict = self.score(materials=m, time=t, efficiency=junction_efficiency)
        return score_dict

    def none_primer_extension_cost(self) -> dict:
        span_to_score_dict = {}
        for span in range(-150, 150):

            span_to_score_dict.setdefault(span, []).append(
                self.score_primer_extension_cost(span)
            )
        return span_to_score_dict

    def single_primer_extension_cost(self) -> dict:
        """Calculate scores for the case where only 1 fragment is extendable by primers

        -------|                fragment 1
            <------             rev Primer
                |--------       fragment
        """
        span_to_score_dict = {}
        for span in range(-150, 150):
            for (ext1, val1) in self.PRIMER_COST_PER_BP.items():
                span_to_score_dict.setdefault(span, []).append(
                    self.score_primer_extension_cost(span, val1[0], ext1)
                )
        return span_to_score_dict

    def dual_primer_extension_cost(self) -> dict:
        """Calculate scores for the case where both fragments are extendable by primers

        -------|                fragment 1
            <------             rev Primer
                 ------>        fwd primer
                    |--------   fragment
        """
        span_to_score_dict = {}
        for span in range(-150, 150):
            for (ext1, val1), (ext2, val2) in combinations(
                self.PRIMER_COST_PER_BP.items(), r=2
            ):
                span_to_score_dict.setdefault(span, []).append(
                    self.score_primer_extension_cost(span, val1[0], ext1, val2[0], ext2)
                )
        return span_to_score_dict

    def synthesis_cost(self, primer_ext_cost_left: dict, primer_ext_cost_right: dict):
        """
        Computes the cost for a large span that supposed to connect two distance fragments.

        1. shifting the synthesized fragment
        2. one, both, or none of the fragments are extendable

        ------------|                 |---------    fragments
                <-------        --------->          primers
                      |-------------|               synthesized fragment

        :return:
        """

        span_to_score_dict = {}

        for span in range(100, 2000, 3):
            span_to_score_dict.setdefault(span, [])
            for gene_size, gene_cost in self.GENE_SYNTHESIS_COST_RANGES.items(3):
                gene_cost = gene_cost[0]
                # e.g. 500bp spanned by 400bp gene would have 100bp span
                #      500 + 0 - 400 == 100
                # e.g. 500bp spanned by 400bp gene with -100 left offset would have 200bp span
                #      500 - -100 - 400 = 200

                max_overlap = gene_size - span
                for left_span in range(0, )
                left_span = -20  # relative to first fragment
                right_span = span - left_span - gene_size
                left_cost = primer_ext_cost_left[left_span]
                right_cost = primer_ext_cost_right[right_span]
                score_dict = self.score(
                    materials=left_cost["materials"]
                    + right_cost["materials"]
                    + gene_cost["base"],
                    time=max(
                        [left_cost["time"], right_cost["time"], gene_cost["time"]]
                    ),
                    efficiency=(
                        left_cost["efficiency"] ** 2 + right_cost["efficiency"] ** 2
                    )
                    ** 0.5,
                )
                span_to_score_dict[span].append(score_dict)
        return span_to_score_dict

    def compute_score(self, materials, time, efficiency):
        if efficiency == 0:
            return self.INF
        m = materials * self.SCORING["materials"]
        t = time * self.SCORING["time"]
        if self.SCORING["efficiency"] == 0:
            e = 1
        else:
            e = 1.0 / efficiency * self.SCORING["efficiency"]
        return (m * t) * e

    def sort_scores(self, d):
        sorted_dict = {}
        for span, score_dict_list in d.items():
            tmp = [dict(_d) for _d in score_dict_list]
            for s in tmp:
                s["score"] = self.compute_score(**s)
            sorted_dict[span] = sorted(tmp, key=lambda x: x["score"])
        return sorted_dict

    def resolve_scores(self, d):
        return {k: v[0] for k, v in self.sort_scores(d).items()}

    # def _optimize_primer(self, primer_choices):
    #     """
    #     Create a dictionary of gap_span to best_extension case.
    #
    #     :param gap_span_start:
    #     :type gap_span_start:
    #     :param gap_span_end:
    #     :type gap_span_end:
    #     :param gap_span_step:
    #     :type gap_span_step:
    #     :param primer_choices:
    #     :type primer_choices:
    #     :return:
    #     :rtype:
    #     """
    #     results = {}
    #
    #     # iterate through gap span
    #     for gap_span in tqdm(
    #         range(*self.JUNCTION_GAP_SPAN_RANGE), desc="optimizing primers"
    #     ):
    #         results[gap_span] = None
    #         for choice in primer_choices:
    #             if not isinstance(choice, tuple):
    #                 choice = (choice,)
    #             extensions = [x["extension"] for x in choice]
    #             costs = [x["cost"] for x in choice]
    #             names = [x["name"] for x in choice]
    #
    #             total_extension = sum(extensions)
    #             total_cost = sum(costs)
    #             junction_efficiency = self._junction_efficiency(
    #                 -gap_span + total_extension
    #             )
    #
    #             # score
    #             score = self.INF
    #             if junction_efficiency > 0:
    #                 score = total_cost / junction_efficiency
    #
    #             _x = results[gap_span]
    #             if _x is None or score < _x["score"]:
    #                 results[gap_span] = {
    #                     "score": score,
    #                     "cost": total_cost,
    #                     "efficiency": junction_efficiency,
    #                     "condition": f"{names}",
    #                     "extension": f"{extensions}",
    #                 }
    #     return results
