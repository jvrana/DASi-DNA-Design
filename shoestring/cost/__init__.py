import seaborn as sns
import pylab as plt
import numpy as np
import pandas as pd


class JunctionCost(object):

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

    day_cost = 20.0
    material_importance = 1.0

    min_anneal = 16
    primers = np.array(
        [
            [16.0, 60.0, 0.0, 0.3, 1.5],
            [45.0, 200.0, 0.0, 0.8, 2.5],
            [min_anneal, min_anneal + 1, 0.0, 0.0, 0.5],  # special case
        ]
    )

    def __init__(self):
        self.min_span = (
            -300
        )  # the minimum spanning distance to evaluate. Means span[0] == -300

        p = []  # cost array
        a = []  # ranges array
        for row in self.primers:
            _a = np.arange(row[0], row[1], dtype=np.int64)
            a.append(_a)
            p.append(np.broadcast_to(row[2:], (_a.shape[0], row[2:].shape[0])))
        a = np.concatenate(a).reshape(-1, 1)
        p = np.concatenate(p)

        assert len(a) == len(p)

        # 2D materials cost matrix
        m = p[:, 1, np.newaxis] * a + p[:, 0, np.newaxis]
        m = m + m.T

        # 2D time cost matrix
        t = p[:, 2, np.newaxis]
        t = np.maximum(t, t.T)

        # the spanning distance
        max_span = 2 * a.max()
        span = np.arange(self.min_span, max_span - self.min_span, 1)

        # extension array
        ext = a - self.min_anneal

        # the span relative to extensions (i.e. the overlap)
        # overlap (sum of extensions - span) for primers of lengths a[x], a[y]
        relative_span = span - (ext + ext.T)[:, :, np.newaxis]
        relative_span = relative_span.swapaxes(2, 0).swapaxes(1, 2)

        # sanity checks
        # span=0, left_size = 200, right_size
        assert relative_span[0].shape == (ext.shape[0], ext.shape[0])
        # span=-300, primer_length=0, primer_length=0
        assert (
            relative_span[0, 0, 0] == self.min_span
        )  # span=-300, primer_length=0, primer_length=0
        assert relative_span[0, 1, 0] == self.min_span - 1
        assert relative_span[0, 0, 1] == self.min_span - 1
        assert relative_span[0, 1, 1] == self.min_span - 2

        # final costs
        e = self.jxn_efficiency[
            np.clip(-relative_span, 0, len(self.jxn_efficiency) - 1)
        ]
        self.xyz_costs = (m * self.material_importance + t * self.day_cost) * 1.0 / e

        self.min_cost_dict = {
            0: self.xyz_costs[:, -1:, -1:].min(axis=(1, 2)),
            1: self.xyz_costs[:, -1:, :-1].min(axis=(1, 2)),
            2: self.xyz_costs[:, :-1, :-1].min(axis=(1, 2)),
        }

    def plot(self):
        df = pd.DataFrame()
        df["span"] = self.min_span + np.arange(len(self.min_cost_dict[0]))
        df["none"] = self.min_cost_dict[0]
        df["one"] = self.min_cost_dict[1]
        df["two"] = self.min_cost_dict[2]
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
        flexibility = []
        span = []
        for i, x in enumerate(self.xyz_costs):

            span.append(i)
            if x.min() != np.Inf:
                opts = np.argwhere(x < x.min() + 10.0)
                flexibility.append(len(opts))
            else:
                flexibility.append(0)

        sns.lineplot(x=np.array(span) + self.min_span, y=flexibility)
        plt.title("Design Flexibility")
        plt.show()

    def junction_cost(self, x, ext=2):
        i = x - self.min_span
        min_cost_per_span = self.min_cost_dict[ext]
        return min_cost_per_span[np.clip(i, 0, len(min_cost_per_span) - 1)]

    def __getitem__(self, span):
        return self.junction_cost(span)


class SynthesisCost(object):
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

    gene_sizes = np.arange(len(gene_synthesis_cost)).reshape(-1, 1)
    gene_costs = gene_synthesis_cost[:, 0].reshape(-1, 1)
    gene_times = gene_synthesis_cost[:, 1].reshape(-1, 1)
    step_size = 10

    def __init__(self, jxn_cost: JunctionCost):
        self.jxn_cost = jxn_cost
        self.cost_dict = None
        self.cost_min_dict = None
        self.span_range = (0, 3000)
        self.left_span_range = (-500, 500)
        self.make_cost_dict()

    def compute_synthesis_costs(self, left_ext=0, right_ext=0):
        sizes = self.gene_sizes[:: self.step_size, :]

        span = np.arange(
            self.span_range[0], self.span_range[1], self.step_size
        ).reshape(1, -1)
        left_span = np.arange(
            self.left_span_range[0], self.left_span_range[1], self.step_size
        ).reshape(-1, 1)[:, np.newaxis]
        left_cost = self.jxn_cost.junction_cost(left_span, ext=left_ext)
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

        right_cost = self.jxn_cost.junction_cost(right_span, ext=right_ext)
        gene_cost = self.gene_costs[sizes]
        gene_time = self.gene_times[sizes]

        # TODO:should just get the material costs here and multiply the efficiencies later
        ext_costs = left_cost + right_cost

        # size, left span, span
        ext_costs = ext_costs.swapaxes(0, 1)
        material = ext_costs + gene_cost
        total = material + gene_time * 20.0
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
        d = {"step": self.step_size}

        d[0] = self.compute_synthesis_costs(0, 0)
        d[1] = self.compute_synthesis_costs(1, 0)
        d[2] = self.compute_synthesis_costs(1, 1)
        self.cost_dict = d

        m = {"step": self.step_size}
        m[0] = d[0].min(axis=(0, 1)).flatten()
        m[1] = d[1].min(axis=(0, 1)).flatten()
        m[2] = d[2].min(axis=(0, 1)).flatten()
        self.cost_min_dict = m

    def get_synthesis_cost(self, span, ext):
        step_size = self.cost_dict["step"]
        return self.cost_min_dict[ext][int(span / step_size)]

    def plot(self):
        df = pd.DataFrame()
        df["e0"] = self.cost_min_dict[0]
        df["e1"] = self.cost_min_dict[1]
        df["e2"] = self.cost_min_dict[2]
        df["span"] = np.arange(self.span_range[0], self.span_range[1], self.step_size)
        df = pd.melt(
            df, id_vars=["span"], value_vars=["e0", "e1", "e2"], value_name="cost"
        )
        print(df)
        sns.lineplot(x="span", y="cost", hue='variable', data=df)
        plt.show()

    def __getitem__(self, span_ext_tuple: tuple):
        return self.get_synthesis_cost(*span_ext_tuple)
