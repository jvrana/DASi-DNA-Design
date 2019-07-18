from shoestring.cost import Cost
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join, abspath, dirname
from shoestring.utils import InfiniteDict

import pytest


@pytest.fixture(scope='function')
def out():
    return join(abspath(dirname(__file__)), "out")


def plot_score_dict(path, score_dict, x, y, prefix=""):
    d = []
    for k, v in score_dict.items():
        _d = dict(v)
        _d["span"] = k
        d.append(_d)
    df = pd.DataFrame(d)
    print(df)
    fig = plt.figure(figsize=(6, 5))
    ax = fig.gca()
    sns.scatterplot(x=x, y=y, data=df, ax=ax)
    # ax.set_ylim(0, 500)
    ax.set_title("{}__{} vs {}".format(prefix, x, y))
    fig.savefig(join(path, "{}__{}vs{}.png".format(prefix, x, y)), format="png")


class TestAndPlotPrimerExtensionCost:
    def test_dual_primer_extension(self, out):
        gac = Cost()
        score_dict = gac.dual_primer_extension_cost()
        score_dict = gac.resolve_scores(score_dict)

        prefix = "DualPrimerExt"
        plot_score_dict(out, score_dict, "span", "materials", prefix)
        plot_score_dict(out, score_dict, "span", "time", prefix)
        plot_score_dict(out, score_dict, "span", "efficiency", prefix)
        plot_score_dict(out, score_dict, "span", "score", prefix)

    def test_single_primer_extension(self, out):
        gac = Cost()
        score_dict = gac.single_primer_extension_cost()
        score_dict = gac.resolve_scores(score_dict)

        prefix = "SinglePrimerExt"
        plot_score_dict(out, score_dict, "span", "materials", prefix)
        plot_score_dict(out, score_dict, "span", "time", prefix)
        plot_score_dict(out, score_dict, "span", "efficiency", prefix)
        plot_score_dict(out, score_dict, "span", "score", prefix)

    def test_none_primer_extension(self, out):
        gac = Cost()
        score_dict = gac.none_primer_extension_cost()
        score_dict = gac.resolve_scores(score_dict)

        prefix = "NoPrimerExt"
        plot_score_dict(out, score_dict, "span", "materials", prefix)
        plot_score_dict(out, score_dict, "span", "time", prefix)
        plot_score_dict(out, score_dict, "span", "efficiency", prefix)
        plot_score_dict(out, score_dict, "span", "score", prefix)


def test_synthesis_cost(out):
    gac = Cost()
    one = gac.resolve_scores(gac.dual_primer_extension_cost())
    one = InfiniteDict(one)
    scores = gac.synthesis_cost(one, one)
    pass
