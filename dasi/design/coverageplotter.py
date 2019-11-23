from typing import Union

import numpy as np
import pylab as plt
from matplotlib.gridspec import GridSpec

from dasi.design import Design
from dasi.design import LibraryDesign
from dasi.utils.sequence_complexity import DNAStats


# TODO: Implement assembly report
#       1. Coverage report for all groups (total and per each)
#       2. Sequence complexity report
#       3. Spanning distance report
#       4. Assembly report
#       5. Tables (assembly table) and reaction table
#       6. Quilt plot for assembly
# TODO: Report print outs (Markdown, HTML, PDF)
class Report:
    def __init__(self, design: Union[Design, LibraryDesign]):
        self.design = design

    @staticmethod
    def plot_coverage(container, query):
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
        gs = GridSpec(3, 2, hspace=0.55)
        ax1 = fig.add_subplot(gs[0, :2])
        ax2 = fig.add_subplot(gs[1, :2], sharex=ax1)
        ax3 = fig.add_subplot(gs[2, :2], sharex=ax1)

        ax1.set_title("Overall Alignment Coverage")
        ax1.set_ylabel("Coverage")
        ax1.set_xlabel("bp")
        ax1.plot(np.sum(data, axis=0))

        ax2.set_title("Alignment Coverage")
        ax2.set_ylabel("Coverage")
        ax2.set_xlabel("bp")
        ax2.plot(data.T)
        ax2.legend(keys, loc="center left", ncol=1, bbox_to_anchor=(1.0, 0.5))

        stats = DNAStats(query + query + query, 14, 20, 20)
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

            y = y[: -(y.shape[0] - x.shape[0])]
            costs_arr.append(y)
            bp_arr.append(x)

        ax3.set_xlabel("bp")
        ax3.set_ylabel("Complexity")
        ax3.set_title("Sequence Complexity".format(window))

        for x, y, l in zip(bp_arr, costs_arr, windows):
            ax3.plot(x, y, label=l)
            ax3.legend(title="window (bp)")

        ax1.set_xlim(0, len(query))
