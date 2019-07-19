# this quickly computes the minimum cost per span for primers
# todo: handle negative span

import numpy as np


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

    primers = np.array([[16.0, 60.0, 0.0, 0.3, 1.5], [45.0, 200.0, 0.0, 0.8, 2.5]])

    def __init__(self):
        self.min_span = -300

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
        ext = a - min(a)

        # overlap (sum of extensions - span) for primers of lengths a[x], a[y]
        relative_span = span - (ext + ext.T)[:, :, np.newaxis]
        relative_span = relative_span.swapaxes(2, 0).swapaxes(1, 2)

        # sanity checks
        assert relative_span[0].shape == (ext.shape[0], ext.shape[0])
        assert relative_span[0, 0, 0] == self.min_span
        assert relative_span[0, 1, 0] == self.min_span - 1
        assert relative_span[0, 0, 1] == self.min_span - 1
        assert relative_span[0, 1, 1] == self.min_span - 2

        # final costs
        e = self.jxn_efficiency[
            np.clip(-relative_span, 0, len(self.jxn_efficiency) - 1)
        ]
        self.xyz_costs = (m * self.material_importance + t * self.day_cost) * 1.0 / e

    def plot(self):
        min_cost_per_span = self.xyz_costs[:, :, :].min(axis=(1, 2))

        import seaborn as sns
        import pylab as plt

        fig = plt.figure(figsize=(6, 5))
        ax = fig.gca()
        sns.scatterplot(
            y=min_cost_per_span,
            x=self.min_span + np.arange(len(min_cost_per_span)),
            ax=ax,
        )
        plt.show()

        flexibility = []
        span = []
        for i, x in enumerate(self.xyz_costs):

            span.append(i)
            if x.min() != np.Inf:
                opts = np.argwhere(x < x.min() + 10.0)
                flexibility.append(len(opts))
            else:
                flexibility.append(0)

        sns.scatterplot(x=np.array(span) + self.min_span, y=flexibility)
        plt.title("Design Flexibility")

    def junction_cost(self, x, ext=2):
        if ext == 2:
            min_cost_per_span = self.xyz_costs[:, :, :].min(axis=(1, 2))
        elif ext == 1:
            min_cost_per_span = self.xyz_costs[:, :1, :].min(axis=(1, 2))
        else:
            min_cost_per_span = self.xyz_costs[:, :1, :1].min(axis=(1, 2))
        i = x - self.min_span
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
        self.make_cost_dict(self.step_size)

    def compute_synthesis_costs(self, left_ext=0, right_ext=0, step_size=5):
        sizes = self.gene_sizes[::step_size, :]

        span = np.arange(0, 3000, step_size).reshape(1, -1)
        left_span = np.arange(-500, 500, step_size).reshape(-1, 1)[:, np.newaxis]
        left_cost = self.jxn_cost.junction_cost(left_span, ext=left_ext)
        right_span = span - sizes - left_span
        right_cost = self.jxn_cost.junction_cost(right_span, ext=right_ext)
        gene_cost = self.gene_costs[sizes]
        gene_time = self.gene_times[sizes]

        # TODO:should just get the material costs here and multiply the efficiencies later
        ext_costs = left_cost + right_cost

        # size, left span, span
        ext_costs = ext_costs.swapaxes(0, 1)
        material = ext_costs + gene_cost
        total = material + gene_time * 20.0
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

    def make_cost_dict(self, step_size):
        d = {"step": step_size}

        d[0] = self.compute_synthesis_costs(0, 0, step_size).min(axis=(0, 1)).flatten()
        d[1] = self.compute_synthesis_costs(1, 0, step_size).min(axis=(0, 1)).flatten()
        d[2] = self.compute_synthesis_costs(1, 1, step_size).min(axis=(0, 1)).flatten()
        self.cost_dict = d
        return d

    def get_synthesis_cost(self, span, ext):
        step_size = self.cost_dict["step"]
        return self.cost_dict[ext][int(span / step_size)]

    def __getitem__(self, span_ext_tuple: tuple):
        return self.get_synthesis_cost(*span_ext_tuple)
