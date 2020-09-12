from typing import List
from typing import TypeVar
from typing import Union

import numpy as np
import pylab as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from dasi.config import Config
from dasi.utils.sequence import DNAStats

DesignType = TypeVar("Design")


# TODO: Implement assembly report
#       1. Coverage report for all groups (total and per each)
#       2. Sequence complexity report
#       3. Spanning distance report
#       4. Assembly report
#       5. Tables (assembly table) and reaction table
#       6. Quilt plot for assembly
# TODO: Report print outs (Markdown, HTML, PDF)


class Report:
    def __init__(self, design: DesignType):
        self.design = design

    def plot_coverage(self, keys: List = None, show: bool = False):
        if not keys:
            keys = self.design.query_keys
        plots = {}
        for qk in keys:
            container = self.design.containers[qk]
            query = self.design.seqdb[qk]
            fig, axes = self.plot_coverage_of_container(container, query)
            plots[qk] = fig, axes
            if show:
                plt.show()
        return plots

    @staticmethod
    def plot_coverage_of_container(container, query):
        arrs = []
        keys = []
        for gkey, groups in container.groups_by_type().items():
            keys.append(gkey)
            a = np.zeros(len(query))
            for group in groups:
                for s in group.query_region.slices():
                    a[s] += 1
            arrs.append(a)
        data = np.vstack(arrs)

        params = {"legend.fontsize": 8}
        plt.rcParams.update(params)

        fig = plt.figure(figsize=(10, 7.5))
        fig.suptitle(query.name, fontsize=16)

        gs = GridSpec(3, 2, hspace=0.55)
        ax1 = fig.add_subplot(gs[0, :2])
        ax2 = fig.add_subplot(gs[1, :2], sharex=ax1)
        ax3 = fig.add_subplot(gs[2, :2], sharex=ax1)

        ax1.set_title("Total Alignment Coverage")
        ax1.set_ylabel("Coverage")
        ax1.set_xlabel("bp")
        ax1.set_yscale("log")
        ax1.plot(np.sum(data, axis=0))

        ax2.set_title("Alignment Coverage")
        ax2.set_ylabel("Coverage")
        ax2.set_xlabel("bp")
        ax2.set_yscale("log")
        ax2.plot(data.T)
        ax2.legend(keys, loc="center left", ncol=1, bbox_to_anchor=(1.0, 0.5))

        stats = DNAStats(str(query.seq) + str(query.seq) + str(query.seq), 14, 20, 20)
        costs_arr = []
        bp_arr = []
        windows = [100, 500, 1000, 2000]
        for window in windows:
            costs = []
            step = min(windows)
            x = np.arange(len(query) - window, len(query) * 2 + window)
            for i in x[::step]:
                costs.append(stats.cost(i, i + window))

            x = x - len(query)
            y = np.repeat(costs, step)

            delta = -(y.shape[0] - x.shape[0])
            if delta == 0:
                delta = None
            y = y[:delta]
            costs_arr.append(y)
            bp_arr.append(x)

        ax3.axhline(
            y=Config.SequenceScoringConfig.complexity_threshold,
            color="k",
            linestyle="-",
        )
        ax3.set_yscale("log")
        ax3.set_xlabel("bp")
        ax3.set_ylabel("Complexity")
        ax3.set_title("Sequence Complexity ({})".format(window))

        for x, y, l in zip(bp_arr, costs_arr, windows):
            ax3.plot(x, y, label=l)
            ax3.legend(title="window (bp)")

        ax1.set_xlim(0, len(query))

        axes = [ax1, ax2, ax3]

        return fig, axes


class AlignmentPlotter:
    def __init__(self, length, starts, ends, types, names, costs, title=""):
        segments = zip(starts, ends, types, names, costs)
        self.segments = list(segments)
        self.length = length

        self.yspacer = 0.3
        self.border = 1
        self.text_spacer = 0.2
        self.title = title
        self.segment_styles = {
            "fragment": {"color": "black", "linewidth": 3},
            "synthesis": {"color": "red", "linewidth": 3},
            "overhang": {"color": "blue", "linewidth": 1},
        }

    def plot(self):

        Y = len([f for f in self.segments if f[2] != "overhang"])

        fig = plt.figure(figsize=(5, self.yspacer * Y))

        ax = fig.gca()
        ax.set_xlim(-self.border, self.length + self.border)
        ax.set_ylim(-self.border, Y + self.border)

        # set border and axis
        ax.axes.get_xaxis().set_visible(True)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_xlabel("bp")
        ax.set_frame_on(False)
        ax.set_title(self.title)
        y = Y
        for i, seg in enumerate(self.segments):
            x1, x2, segtype, name, cost = seg
            if "overlap" in name or "overlap" in segtype:
                segtype = "overhang"
            elif "synthesis" in name or "synthesis" in segtype:
                segtype = "synthesis"
            else:
                segtype = "fragment"
            cost = "{:.1f}".format(float(cost))
            name = "{} - {}".format(Y - y, name)
            segstyle = self.segment_styles[segtype]
            if x1 < x2:
                line = Line2D([seg[0], seg[1]], [y, y], zorder=1, **segstyle)

                # add cpst
                xtext = (seg[0] + seg[1]) / 2.0
                ax.text(
                    xtext,
                    y + self.text_spacer,
                    str(cost),
                    horizontalalignment="center",
                    verticalalignment="bottom",
                )

                # add label
                ax.text(
                    self.length,
                    y,
                    name,
                    horizontalalignment="left",
                    verticalalignment="center",
                )
                ax.add_line(line)
                y -= 1
            else:
                if segtype == "overhang":
                    if i == len(self.segments) - 1:
                        y2 = Y
                    else:
                        y2 = y + 1
                    line1 = Line2D([x1, x2], [y, y2], zorder=0, **segstyle)
                    line2 = Line2D([x2, x1], [y, y2], zorder=0, **segstyle)

                    xtext = (x1 + x2) / 2.0
                    ax.text(
                        xtext,
                        y - 2 * self.text_spacer,
                        str(cost),
                        horizontalalignment="center",
                        verticalalignment="top",
                    )

                    ax.add_line(line1)
                    ax.add_line(line2)
                else:
                    line1 = Line2D([x1, self.length], [y, y], zorder=1, **segstyle)
                    line2 = Line2D([0, x2], [y, y], zorder=1, **segstyle)

                    # add cost
                    if self.length - x1 > x2:
                        xtext = (self.length - x1) / 2.0
                    else:
                        xtext = x2 / 2.0
                    ax.text(
                        xtext,
                        y + self.text_spacer,
                        str(cost),
                        horizontalalignment="center",
                        verticalalignment="bottom",
                    )

                    # add label
                    ax.text(
                        self.length,
                        y,
                        name,
                        horizontalalignment="left",
                        verticalalignment="center",
                    )

                    ax.add_line(line1)
                    ax.add_line(line2)
        return fig, ax
