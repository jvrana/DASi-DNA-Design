import seaborn as sns
import pylab as plt
import numpy as np
import pandas as pd

# TODO: backtrack options


class JxnParams(object):

    # jxn_efficiency[30:40] = 0.9" means for 30 <= bp < 40, junction efficiency is 90%
    jxn_efficiency = np.zeros(301, dtype=np.float)
    jxn_efficiency[:10] = 0.0
    jxn_efficiency[10:20] = 0.1
    jxn_efficiency[20:30] = 0.8
    jxn_efficiency[30:40] = 0.9
    jxn_efficiency[40:50] = 0.8
    jxn_efficiency[50:100] = 0.75
    jxn_efficiency[100:120] = 0.5
    jxn_efficiency[120:150] = 0.3
    jxn_efficiency[150:250] = 0.1
    jxn_efficiency[250:300] = 0.0

    min_jxn_span = -300     # the minimum spanning junction we evaluate
    min_anneal = 16         # the minimum annealing of each primer
    primers = np.array(
        [
            [16.0, 60.0, 0.0, 0.3, 1.5],
            [45.0, 200.0, 0.0, 0.8, 2.5],
            [min_anneal, min_anneal + 1, 0.0, 0.0, 0.5],  # special case
        ]
    )
    primer_rows = ["IDTPrimer", "IDTUltramer", "NoPrimer(Free)"]
    primer_cols = ["min_bp", "max_bp", "base_cost", "bp_cost", "time (days)"]

    # sanity check
    assert len(primer_rows) == len(primers)
    assert len(primer_cols) == len(primers[0])


class Slicer(object):

    def __getitem__(self, item):
        return item
    
slicer = Slicer()


class SynParams(object):
    gene_synthesis_cost = np.zeros((2250, 2))
    d = {
        (0, 1): {"base": 0.0, "time": 0},
        (1, 100): {"base": np.Inf, "time": np.Inf},
        (100, 500): {"base": 89.0, "time": 3.0},
        (500, 750): {"base": 129.0, "time": 3.0},
        (750, 1000): {"base": 149.0, "time": 4.0},
        (1000, 1250): {"base": 209.0, "time": 7.0},
        (1250, 1500): {"base": 249.0, "time": 7.0},
        (1500, 1750): {"base": 289.0, "time": 7.0},
        (1750, 2000): {"base": 329.0, "time": 7.0},
        (2000, 2250): {"base": 399.0, "time": 7.0},
    }
    for k, v in d.items():
        gene_synthesis_cost[k[0] : k[1]] = np.array([v["base"], v["time"]])

    synthesis_span_range = (0, 3000)
    gene_sizes = np.arange(len(gene_synthesis_cost)).reshape(-1, 1)
    gene_costs = gene_synthesis_cost[:, 0].reshape(-1, 1)
    gene_times = gene_synthesis_cost[:, 1].reshape(-1, 1)
    synthesis_step_size = 10
    synthesis_left_span_range = (-500, 500)


class CostParams(object):

    time = 50.0  # dollar cost of 24 hours
    material = 1.0  # multiply material cost by this amount


class JunctionCost(object):
    def __init__(self):
        min_span = JxnParams.min_jxn_span
        max_span = JxnParams.primers[:, 1].max() * 2 - min_span
        self.span = np.arange(min_span, max_span, dtype=np.int64)

        p = []  # cost array
        a = []  # ranges array
        for row in JxnParams.primers:
            _a = np.arange(row[0], row[1], dtype=np.int64)
            a.append(_a)
            p.append(np.broadcast_to(row[2:], (_a.shape[0], row[2:].shape[0])))
        a = np.concatenate(a).reshape(-1, 1)
        self.primer_lengths_array = a
        p = np.concatenate(p)

        assert len(a) == len(p)

        # 2D materials cost matrix
        m = p[:, 1, np.newaxis] * a + p[:, 0, np.newaxis]
        m = m + m.T

        # 2D time cost matrix
        t = p[:, 2, np.newaxis]
        t = np.maximum(t, t.T)

        # extension array
        ext = a - JxnParams.min_anneal

        # the span relative to extensions (i.e. the overlap)
        # overlap (sum of extensions - span) for primers of lengths a[x], a[y]
        relative_span = self.span - (ext + ext.T)[:, :, np.newaxis]
        relative_span = relative_span.swapaxes(2, 0).swapaxes(1, 2)

        # sanity checks
        # span=0, left_size = 200, right_size
        assert relative_span[0].shape == (ext.shape[0], ext.shape[0])
        # span=-300, primer_length=0, primer_length=0
        assert (
            relative_span[0, 0, 0] == self.span.min()
        )  # span=-300, primer_length=0, primer_length=0
        assert relative_span[0, 1, 0] == self.span.min() - 1
        assert relative_span[0, 0, 1] == self.span.min() - 1
        assert relative_span[0, 1, 1] == self.span.min() - 2

        # final costs
        e = JxnParams.jxn_efficiency[
            np.clip(-relative_span, 0, len(JxnParams.jxn_efficiency) - 1)
        ]
        self.xyz_labels = ["span", "left_ext", "right_ext"]
        self.cost_matrix = (m * CostParams.material + t * CostParams.time) * 1.0 / e

        self.slice_dict = {
            (0, 0): slicer[:, -1:, -1:],
            (0, 1): slicer[:, -1:, :-1],
            (1, 0): slicer[:, :-1, -1:],
            (1, 1): slicer[:, :-1, :-1],
        }
        self.cost_dict = {k: self.cost_matrix[self.slice_dict[k]] for k in self.slice_dict}
        self.min_cost_dict = {k: v.min(axis=(1, 2)) for k, v in self.cost_dict.items()}

    def plot(self):
        df = pd.DataFrame()
        df["span"] = self.span
        df["none"] = self.min_cost_dict[(0, 0)]
        df["one"] = self.min_cost_dict[(1, 0)]
        df["two"] = self.min_cost_dict[(1, 1)]
        df = pd.melt(
            df, id_vars=["span"], value_vars=["none", "one", "two"], value_name="cost"
        )

        print(df.columns)
        fig = plt.figure(figsize=(6, 5))
        ax = fig.gca()
        sns.lineplot(y="cost", x="span", hue="variable", data=df, ax=ax)
        ax.set_title("cost vs spanning distance (bp)")
        plt.show()

    def plot_design_flexibility(self):
        """Makes a plot of the design flexibility for a given bp span"""
        options_arr = []
        for x in self.cost_matrix:
            if x.min() != np.Inf:
                opts = np.argwhere(x < x.min() + 10.0)
                options_arr.append(len(opts))
            else:
                options_arr.append(0)

        df = pd.DataFrame()
        df["span"] = self.span
        df["options"] = options_arr

        fig = plt.figure(figsize=(6, 5))
        ax = fig.gca()

        sns.lineplot(x="span", y="options", ax=ax, data=df)
        plt.title("Design Flexibility")
        plt.show()

    # def options(self, span, ext, max_cost):
    #     m = self.cost_dict[ext][span]
    #     print(m.min())
    #     indices = np.argwhere(m <= max_cost)
    #     return indices
    #
    # def _span_to_index(self, span):
    #     i = span - self.span.min()
    #     i = np.clip(i, 0, self.cost_dict.shape[0] - 1)
    #     return i
    #
    #
    #
    # def cost_matrix(self, span, ext):
    #     return self.cost_dict[ext][self._span_to_index(span)]

    def cost(self, span, ext):
        i = span - self.span.min()
        i = np.clip(i, 0, self.min_cost_dict[ext].shape[0] - 1)
        return self.min_cost_dict[ext][i]


class SynthesisCost(object):
    def __init__(self, jxn_cost: JunctionCost):
        self.jxn_cost = jxn_cost
        self.cost_dict = None
        self.cost_min_dict = None
        self.span = np.arange(
            SynParams.synthesis_span_range[0],
            SynParams.synthesis_span_range[1],
            SynParams.synthesis_step_size,
        )
        self.sizes = SynParams.gene_sizes[:: SynParams.synthesis_step_size, :]
        self.make_cost_dict()

    def compute_synthesis_costs(self, left_ext=0, right_ext=0):
        sizes = self.sizes
        span = self.span.reshape(1, -1)

        left_span = np.arange(
            SynParams.synthesis_left_span_range[0],
            SynParams.synthesis_left_span_range[1],
            SynParams.synthesis_step_size,
        ).reshape(-1, 1)[:, np.newaxis]
        left_cost = self.jxn_cost.cost(left_span, ext=(left_ext, 0))
        right_span = span - sizes - left_span

        # sanity check
        # left_span[0] == -500
        # gene_sizes[0] * step_size == 0
        # span[0] * step_size == 0
        assert right_span[0, 100, 100] == 500

        # left_span[10] == -400
        # gene_sizes[100] == 1000
        # span[50] == 500
        assert right_span[10, 100, 50] == -100

        right_cost = self.jxn_cost.cost(right_span, ext=(0, right_ext))
        gene_cost = SynParams.gene_costs[sizes]
        gene_time = SynParams.gene_times[sizes]

        # TODO:should just get the material costs here and multiply the efficiencies later
        ext_costs = left_cost + right_cost

        # size, left span, span
        ext_costs = ext_costs.swapaxes(0, 1)
        material = ext_costs + gene_cost
        total = material * CostParams.material + gene_time * CostParams.time
        # axes are [gene_size, left_span, span]
        return total

    def optimize_step_size(self, s, ds):
        step_size = s
        delta_step = ds
        y1 = self.compute_synthesis_costs(0, 0, step_size).min(axis=(0, 1)).flatten()
        x1 = np.arange(len(y1)) * step_size

        y2 = (
            self.compute_synthesis_costs(0, 0, step_size - delta_step)
            .min(axis=(0, 1))
            .flatten()
        )
        x2 = np.arange(len(y2)) * (step_size - delta_step)
        y2_interp = np.interp(x1, x2, y2)
        diff = ((y1 - y2_interp) ** 2).sum() / len(y1)
        return diff

    def make_cost_dict(self):
        d = {}
        d[(0, 0)] = self.compute_synthesis_costs(0, 0)
        d[(1, 0)] = self.compute_synthesis_costs(1, 0)
        d[(0, 1)] = self.compute_synthesis_costs(0, 1)
        d[(1, 1)] = self.compute_synthesis_costs(1, 1)
        self.cost_dict = d
        self.cost_min_dict = {
            k: v.min(axis=(0, 1)).flatten() for k, v in self.cost_dict.items()
        }
        self.cost_dict["step"] = SynParams.synthesis_step_size
        self.cost_min_dict["step"] = SynParams.synthesis_step_size

    def cost(self, span, ext):
        assert self.span[0] == 0
        step_size = self.cost_dict["step"]
        i = np.array(span / step_size, dtype=np.int64)
        i = np.clip(i, 0, len(self.span) - 1)
        return self.cost_min_dict[ext][i]

    def plot(self):
        df = pd.DataFrame()
        df["e0"] = self.cost_min_dict[(0, 0)]
        df["e1"] = self.cost_min_dict[(1, 0)]
        df["e2"] = self.cost_min_dict[(1, 1)]
        df["span"] = np.arange(
            SynParams.synthesis_span_range[0],
            SynParams.synthesis_span_range[1],
            SynParams.synthesis_step_size,
        )
        df = pd.melt(
            df, id_vars=["span"], value_vars=["e0", "e1", "e2"], value_name="cost"
        )
        fig = plt.figure(figsize=(6, 5))
        ax = fig.gca()
        ax.set_ylim(0, 1000)
        sns.lineplot(x="span", y="cost", hue="variable", data=df, ax=ax)
        plt.show()


class SpanCost(object):
    def __init__(self):
        self.junction_cost = JunctionCost()
        self.synthesis_cost = SynthesisCost(self.junction_cost)

    def cost(self, span: int, ext: int) -> float:
        x = np.array(
            [self.junction_cost.cost(span, ext), self.synthesis_cost.cost(span, ext)]
        )
        return x.min(axis=0)
