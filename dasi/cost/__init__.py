"""Cost.

.. module:: dasi.cost

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    params
    utils
"""

from abc import ABC, abstractmethod
from functools import partial

import numpy as np
import pandas as pd
import seaborn as sns

from dasi.utils import NumpyDataFrame
from .params import PrimerParams, SynthesisParams, Globals
from .utils import lexargmin, slicer, df_to_np_ranged
import msgpack
from typing import Union, Tuple

slice_dict = {
    (0, 0): slicer[:, :1, :1],
    (0, 1): slicer[:, :1, 1:],
    (1, 0): slicer[:, 1:, :1],
    (1, 1): slicer[:, 1:, 1:],
}


def encoder(obj):
    """msgpack encoder for cost functions.

    :param obj:
    :return:
    """
    if isinstance(obj, NumpyDataFrame):
        return {"__numpydataframe__": True, "data": obj.data}
    elif isinstance(obj, SpanCost):
        return {
            "__spancost__": True,
            "cost_dict": {k: v for k, v in obj.cost_dict.items()},
            "span": obj.span,
        }
    return obj


def decoder(obj):
    """msgpack decoder for cost functions.

    :param obj:
    :return:
    """
    if b"__numpydataframe__" in obj:
        data = obj[b"data"]
        data = {k.decode(): v for k, v in data.items()}
        obj = NumpyDataFrame(data=data)
    elif b"__spancost__" in obj:
        cost_dict = {tuple(k): v for k, v in obj[b"cost_dict"].items()}
        span = obj[b"span"]
        obj = SpanCost.__new__(SpanCost)
        obj.cost_dict = cost_dict
        obj.span = span
    return obj


class CostBuilder(ABC):
    def __init__(self, span):
        self.cost_dict = {}
        self.span = span

    @classmethod
    @abstractmethod
    def from_params(cls, *args):
        pass

    @abstractmethod
    def compute(self):
        """Before the numpy-based cost calculations."""
        pass

    def to_df(self):
        """Convert the cost builder to a pandas.DataFrame."""
        dfs = []
        for ext, jxn in self.cost_dict.items():
            df = jxn.to_df()
            df["condition"] = [ext] * len(df)
            dfs.append(df)
        df = pd.concat(dfs)
        return df

    @property
    def plot(self):
        """Plot attributes across 'span' from the pandas.DataFrame."""
        return partial(
            sns.lineplot, data=self.to_df(), x="span", y="cost", hue="condition"
        )

    # TODO: what if there is a gap in the span?
    def cost(
        self, bp: Union[np.ndarray, int], ext: Tuple[bool, bool], invalidate=False
    ):
        """

        :param bp:
        :param ext:
        :param invalidate:
        :return:
        """
        if isinstance(bp, int):
            bp = np.array([bp])
        _span = self.span.flatten()
        index = np.array(bp - _span.min(), dtype=np.int64)

        # clipped index
        clipped = np.clip(index, 0, len(_span) - 1)

        # get junction
        jxn = self.cost_dict[ext][clipped]

        if invalidate:
            # invalidated indices
            invalid_i1 = index < 0
            invalid_i2 = index >= len(_span)

            # invalidate extremes
            for invalid in [invalid_i1, invalid_i2]:
                jxn.data["cost"][invalid] = np.inf
                jxn.data["efficiency"][invalid] = 0.0

        # add span
        jxn.col["span"] = bp
        return jxn

    def __call__(self, bp, ext):
        return self.cost(bp, ext)


class PrimerCostBuilder(CostBuilder):
    def __init__(
        self,
        pdf: pd.DataFrame,
        edf: pd.DataFrame,
        min_anneal: int,
        time_cost: float,
        material_mod: float,
        min_span: int,
    ):
        self.primer_df = pdf
        self.eff_df = edf
        self.time_cost = time_cost
        self.min_anneal = min_anneal
        max_span = self.primer_df["max"].max() * 2 - min_span
        span = np.arange(min_span, max_span, dtype=np.int32)
        super().__init__(span)
        self.material_modifier = material_mod
        self.compute()

    @classmethod
    def from_params(cls, primer_params: PrimerParams):
        return cls(
            pdf=primer_params.primer_df,
            edf=primer_params.eff_df,
            min_anneal=primer_params.min_anneal,
            time_cost=primer_params.time_cost,
            material_mod=primer_params.material_modifier,
            min_span=primer_params.min_span,
        )

    def compute(self):
        span = self.span

        # span, base cost, cost per bp, time (days)
        p = df_to_np_ranged(
            "min",
            "max",
            self.primer_df,
            cols=["base cost", "cost per bp", "time (days)"],
            dtype=np.float64,
        )

        # flattened extension array
        ext = p[:, 0].reshape(-1, 1) - self.min_anneal
        ext = ext.astype(np.int32)

        # relative span (i.e. the overlap)
        rel_span = span[:, np.newaxis, np.newaxis] - (ext + ext.T)[np.newaxis, :, :]

        # efficiency, the same shape as rel_span
        eff_arr = df_to_np_ranged("min", "max", self.eff_df, dtype=np.float64)[:, 1]
        eff = eff_arr[np.clip(-rel_span, 0, len(eff_arr) - 1)]

        # material cost
        m = p[:, 0, np.newaxis] * p[:, 2, np.newaxis] + p[:, 1, np.newaxis]
        t = p[:, 3, np.newaxis]
        t = np.maximum(t, t.T)
        x = m * self.material_modifier + t * self.time_cost
        material_cost = x + x.T

        # cost
        cost = material_cost / eff
        cost[np.where(np.isnan(cost))] = np.inf

        for slice_index, slice_obj in slice_dict.items():
            s_eff = eff[slice_obj]
            s_cost = cost[slice_obj]
            s_mat = material_cost[slice_obj[1], slice_obj[2]]

            idx = lexargmin((s_eff, s_cost), axis=0)
            self.cost_dict[slice_index] = NumpyDataFrame(
                dict(
                    span=span[idx[0]],
                    cost=s_cost[idx],
                    efficiency=s_eff[idx],
                    material=s_mat[idx[1], idx[2]],
                    left_ext=ext[idx[1]],
                    right_ext=ext[idx[2]],
                    time=t[idx[1], idx[2]],
                ),
                apply=np.squeeze,
            )


class SynthesisCostBuilder(CostBuilder):
    def __init__(
        self,
        sdf: pd.DataFrame,
        primer_cost: PrimerCostBuilder,
        time_cost: float,
        material_modifier: float,
        step_size=10,
        left_span_range=(-500, 500),
    ):
        self.synthesis_df = sdf
        self.step_size = step_size
        self.material_modifier = (material_modifier,)
        self.lspanrange = left_span_range
        self.primer_cost = primer_cost
        self.time_cost = time_cost
        super().__init__(None)
        # self.logger = logger(self)
        self.compute()

    @classmethod
    def from_params(cls, syn_params: SynthesisParams, primer_cost: PrimerCostBuilder):
        return cls(
            sdf=syn_params.synthesis_df,
            primer_cost=primer_cost,
            time_cost=syn_params.time_cost,
            material_modifier=syn_params.material_modifier,
            step_size=syn_params.step_size,
            left_span_range=syn_params.left_span_range,
        )

    def compute(self):
        syn = df_to_np_ranged("min", "max", self.synthesis_df, dtype=np.float32)
        gene_sizes = syn[:, 0].reshape(-1, 1)[:: self.step_size, :].astype(np.int32)
        gene_costs = syn[:, 1][gene_sizes]
        gene_times = syn[:, 2][gene_sizes]

        self.span = syn[:, 0].reshape(1, -1).astype(np.int32)

        left_span = np.arange(self.lspanrange[0], self.lspanrange[1], self.step_size)[
            ..., np.newaxis, np.newaxis
        ]

        for i, j in [(0, 0), (1, 0), (0, 1), (1, 1)]:
            jxn = self._compute(gene_costs, gene_sizes, gene_times, i, j, left_span)
            self.cost_dict[(i, j)] = jxn

    def _compute(self, gene_costs, gene_sizes, gene_times, i, j, left_span):
        # extension conditions
        left_ext = (i, 0)
        right_ext = (0, j)
        # left primer
        left_jxn = self.primer_cost(left_span, ext=left_ext)
        left_eff = left_jxn.data["efficiency"]
        left_material = left_jxn.data["material"]

        # right primer
        right_span = self.span - gene_sizes - left_span
        right_jxn = self.primer_cost(right_span, ext=right_ext)
        right_eff = right_jxn.data["efficiency"]
        right_material = right_jxn.data["material"]
        ext_material = left_material + right_material
        ext_eff = np.multiply(left_eff, right_eff)

        # swap axes
        # span, size, left_span
        ext_material = ext_material.swapaxes(0, 2)
        ext_eff = ext_eff.swapaxes(0, 2)
        syn_eff = ext_eff * 1.0  # here place probability of success for gene synthesis
        # could even use sequence to compute this later???
        syn_material_cost = (
            ext_material + gene_costs[np.newaxis, ...] * self.material_modifier
        )
        syn_time_cost = gene_times * self.time_cost
        syn_total_cost = (syn_material_cost + syn_time_cost[np.newaxis, ...]) / syn_eff
        idx = lexargmin((syn_eff, syn_total_cost), axis=0)

        _gcosts = gene_costs[idx[1]]
        _span = np.squeeze(self.span)[idx[0]]

        gene_df = NumpyDataFrame(
            dict(
                cost=_gcosts,
                material=_gcosts,
                efficiency=np.ones(idx[0].shape[0]),
                size=gene_sizes[idx[1]],
            ),
            apply=np.squeeze,
        )

        gap_df = NumpyDataFrame(
            dict(
                span=_span,
                cost=syn_total_cost[idx],
                efficiency=syn_eff[idx],
                material=syn_material_cost[idx],
                lshift=left_span[idx[2]],
            ),
            apply=np.squeeze,
        )

        gap_df.update(left_jxn[idx[2]].apply(np.squeeze).prefix("lprimer_"))
        gap_df.update(right_jxn[idx[2], idx[1], idx[0]].prefix("rprimer_"))
        gap_df.update(gene_df.prefix("gene_"))

        return gap_df

    def cost(self, bp, ext, invalidate=True):
        return super().cost(bp, ext, invalidate=invalidate)

    def __call__(self, bp, ext):
        return self.cost(bp, ext)


class SpanCost(CostBuilder):
    def __init__(self, syn_cost):
        self.syn_cost = syn_cost
        self.primer_cost = self.syn_cost.primer_cost
        x = [
            syn_cost.span.min(),
            syn_cost.span.max(),
            self.primer_cost.span.min(),
            self.primer_cost.span.max(),
        ]
        span = np.arange(min(x), max(x))
        super().__init__(span)
        self.compute()

    @classmethod
    def from_params(cls):
        pass

    @classmethod
    def default(cls):
        primer_cost = PrimerCostBuilder.from_params(PrimerParams)
        syn_cost = SynthesisCostBuilder.from_params(SynthesisParams, primer_cost)
        return cls(syn_cost)

    def compute(self):
        def choose(a, i):
            return np.choose(i, a)

        for s in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            # numpy data frames for primer cost and syn cost over span
            df1 = self.primer_cost(self.span, s)
            df2 = self.syn_cost(self.span, s)

            # determine the indices of the min cost (0=primer, 1=syn)
            c1 = df1.data["cost"]
            c2 = df2.data["cost"]
            c3 = np.stack((c1, c2), axis=1)
            y = c3.argmin(axis=1)

            # select between primer_cost and syn_cost based on the min cost
            df4 = NumpyDataFrame.group_apply(
                (df1, df2), choose, i=y, _fill_value=np.nan
            )
            self.cost_dict[s] = df4

    def cost(self, bp, ext, invalidate=True):
        return super().cost(bp, ext, invalidate=invalidate)

    def dumpb(self) -> bytes:
        return msgpack.packb(self, default=encoder, use_bin_type=True)

    @classmethod
    def loadb(cls, s: bytes):
        return msgpack.unpackb(s, object_hook=decoder, raw=True, use_list=False)

    def dump(self, path: str):
        with open(path, "wb") as f:
            f.write(self.dumpb())

    @classmethod
    def load(cls, path: str):
        with open(path, "rb") as f:
            return cls.loadb(f.read())
