import json
from abc import ABC
from abc import abstractmethod
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Tuple
from typing import Union

import msgpack
import numpy as np
import pandas as pd

from .utils import df_to_np_ranged
from .utils import lexargmin
from .utils import slicer
from dasi.exceptions import DasiCostParameterValidationError
from dasi.schemas import Schemas
from dasi.schemas import validate_with_schema
from dasi.utils import NumpyDataFrame

here = abspath(dirname(__file__))


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
    elif isinstance(obj, PrimerCostModel):
        return (
            {
                "__primercostmodel__": True,
                "cost_dict": {k: v for k, v in obj.cost_dict.items()},
                "span": obj.span,
            },
        )
    elif isinstance(obj, SynthesisCostModel):
        return {
            "__synthesiscostmodel__": True,
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
    elif b"__primercostmodel__" in obj:
        cost_dict = {tuple(k): v for k, v in obj[b"cost_dict"].items()}
        span = obj[b"span"]
        obj = PrimerCostModel.__new__(PrimerCostModel)
        obj.cost_dict = cost_dict
        obj.span = span
    elif b"__synthesiscostmodel__" in obj:
        cost_dict = {tuple(k): v for k, v in obj[b"cost_dict"].items()}
        span = obj[b"span"]
        obj = SynthesisCostModel.__new__(SynthesisCostModel)
        obj.cost_dict = cost_dict
        obj.span = span
    elif b"__spancost__" in obj:
        cost_dict = {tuple(k): v for k, v in obj[b"cost_dict"].items()}
        span = obj[b"span"]
        obj = SpanCost.__new__(SpanCost)
        obj.cost_dict = cost_dict
        obj.span = span
    return obj


def _add_special_primer_case(params):
    primer_cost_data = params["primer_cost"]["data"]
    primer_cost_columns = params["primer_cost"]["columns"]
    i = primer_cost_columns.index("min")
    j = primer_cost_columns.index("max")
    k = primer_cost_columns.index("name")

    new = [0] * len(primer_cost_columns)

    new[k] = "Non-primer (special case)"
    new[i] = params["primer_min_anneal"]
    new[j] = params["primer_min_anneal"] + 1
    if new not in primer_cost_data:
        primer_cost_data.insert(0, new)


def _replace_inf(params):
    if isinstance(params, list):
        for i, p in enumerate(params):
            if p == "inf":
                params[i] = np.inf
        for p in params:
            _replace_inf(p)
    elif isinstance(params, dict):

        for k, v in params.items():
            if v == "inf":
                params[k] = np.inf
            else:
                _replace_inf(v)


def validate_params(params):
    schema = Schemas.cost_parameters_schema
    validate_with_schema(params, schema, reraise_as=DasiCostParameterValidationError)


default_parameters_path = join(here, "default_parameters.json")


def open_params(path=None):
    if path is None:
        path = default_parameters_path
    with open(path) as f:
        params = json.load(f)
    return params


class CostBuilder(ABC):
    """Abstract base class for building a cost model."""

    def __init__(self, span: np.ndarray):
        """Initialize cost builder from a 'span', or a list of indices for
        which the model is valid.

        :param span: the list of indices (span) for which this cost model is evaluated.
        """
        self.cost_dict = {}
        self.span = span

    @staticmethod
    def _fix_params(params):
        validate_params(params)
        _add_special_primer_case(params)
        _replace_inf(params)
        return params

    @classmethod
    @abstractmethod
    def from_json(cls, params: dict, override_span: np.ndarray = None):
        """Initialize from parameters."""
        cls._fix_params(params)

    @classmethod
    def open(cls, path: str = None):
        params = open_params(path)
        return cls.from_json(params)

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

    # @property
    # def plot(self):
    #     """Plot attributes across 'span' from the pandas.DataFrame."""
    #     return partial(
    #         sns.lineplot, data=self.to_df(), x="span", y="cost", hue="condition"
    #     )

    # TODO: what if there is a gap in the span?
    def cost(
        self, bp: Union[np.ndarray, int], ext: Tuple[int, int], invalidate=False
    ) -> NumpyDataFrame:
        """Return the :class:`NumpyDataFrame <dasi.usil.NumpyDataFrame>` across
        the provided span.

        :param bp: either a integer or np.ndarray of integers to compute.
        :param ext: The extension design parameters (bool, bool) that represent
                    whether the left or right primer is 'extendable'. Primers
                    that have been provided are not extendable. If we had,
                    for example, an existing right primer and flexibility to
                    design the left primer, this would be `(1, 0)`.
        :param invalidate: Whether to invalidate indices that go beyond the provided
                            span for the cost builder (`CostModelBase.span`). Any
                            invalid indices cost will be set to `np.inf` and efficiency
                            to `0.0`
        :return: NumpyDataFrame representing the cost of from the input span.
        """
        if isinstance(bp, int):
            bp = np.array([bp])
        _span = self.span.flatten()
        index = np.array(bp - _span.min(), dtype=np.int64)

        # clipped index
        clipped = np.clip(index, 0, len(_span) - 1)
        if isinstance(clipped, int):
            index = np.array([clipped])

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
        jxn.col["span"] = np.array(bp)
        return jxn

    def __call__(
        self, bp: Union[int, np.ndarray], ext: Tuple[int, int]
    ) -> NumpyDataFrame:
        """Short hand for calling :meth:`cost`"""
        return self.cost(bp, ext)


class PrimerCostModel(CostBuilder):
    def __init__(
        self,
        pdf: pd.DataFrame,
        edf: pd.DataFrame,
        min_anneal: int,
        time_cost: float,
        material_mod: float,
        min_span: int,
        input_data: dict = None,
        override_span: np.ndarray = None,
    ):
        """Initialize the cost model for primer designs.

        :param pdf: a pandas DataFrame representing costs for primers.
                    Example parameters found at
                    :attr:`dasi.cost.params.PrimerParams.primer_df`
        :param edf: a pandas DataFrame representing efficiencies for spans.
                    Example parameters found at
                    :attr:`dasi.cost.params.PrimerParams.eff_df`
        :param min_anneal: the number of bases to consider for annealing a primer.
        :param time_cost: the cost to produce a primer
        :param material_mod:
        :param min_span:
        """
        self.input_data = input_data
        self.primer_df = pdf
        self.eff_df = edf
        self.time_cost = time_cost
        self.min_anneal = min_anneal
        max_span = self.primer_df["max"].max() * 2 - min_span
        if override_span is not None:
            span = override_span
        else:
            span = np.arange(min_span, max_span, dtype=np.int32)
        super().__init__(span)
        self.material_modifier = material_mod
        self.compute()

    @classmethod
    def from_json(cls, params, override_span: np.ndarray = None) -> "PrimerCostModel":
        """Load from :class:`dasi.cost.params.PrimerParams`.

        :param primer_params: parameters
        :return: PrimerCostModel
        """
        super().from_json(params)
        return PrimerCostModel(
            pdf=pd.DataFrame(
                params["primer_cost"]["data"], columns=params["primer_cost"]["columns"]
            ),
            edf=pd.DataFrame(
                params["primer_efficiency"]["data"],
                columns=params["primer_efficiency"]["columns"],
            ),
            time_cost=params["global_time_cost"],
            material_mod=params["global_material_modifier"],
            min_span=params["_primer_min_span"],
            min_anneal=params["primer_min_anneal"],
            input_data=params,
            override_span=override_span,
        )

    # TODO: IMPORTANT!! For extreme values, for extensions material cost should NOT
    #       BE ZERO
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

        slice_dict = {
            (0, 0): slicer[:, :1, :1],
            (0, 1): slicer[:, :1, 1:],
            (1, 0): slicer[:, 1:, :1],
            (1, 1): slicer[:, 1:, 1:],
        }

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


class SynthesisCostModel(CostBuilder):
    def __init__(
        self,
        sdf: pd.DataFrame,
        primer_cost: PrimerCostModel,
        time_cost: float,
        material_modifier: float,
        step_size=10,
        left_span_range=(-500, 500),
        input_data: dict = None,
        override_span: np.ndarray = None,
    ):
        """

        :param sdf: synthesis parameter pandas dataframe. Example parameters found
                    at :attr:`dasi.cost.params.SynthesisParams.synthesis_df`
        :param primer_cost: primer cost builder
        :param time_cost: time cost per day
        :param material_modifier: material multiplier (default: 1.0)
        :param step_size: step size to consider (default: 10)
        :param left_span_range: span range for left primer to consider
                (default: -500, 500)
        :param input_data: optional input data used to created the cost model
        :param override_span: optionally, only compute over the provided span (for
                debugging/testing)
        """
        self.input_data = input_data
        self.synthesis_df = sdf
        self.step_size = step_size
        self.material_modifier = (material_modifier,)
        self.lspanrange = left_span_range
        self.primer_cost = primer_cost
        self.time_cost = time_cost
        self.syn = df_to_np_ranged("min", "max", self.synthesis_df, dtype=np.float32)
        if override_span is not None:
            span = override_span
        else:
            span = self.syn[:, 0].reshape(1, -1).astype(np.int32)
        super().__init__(span)
        self.compute()

    @classmethod
    def from_json(
        cls,
        params: dict,
        primer_cost_model: PrimerCostModel,
        override_span: np.ndarray = None,
    ):
        super().from_json(params)
        return SynthesisCostModel(
            sdf=pd.DataFrame(
                params["synthesis_cost"]["data"],
                columns=params["synthesis_cost"]["columns"],
            ),
            primer_cost=primer_cost_model,
            time_cost=params["global_time_cost"],
            material_modifier=params["global_material_modifier"],
            step_size=params["_synthesis_step_size"],
            left_span_range=params["_synthesis_left_span_range"],
            override_span=override_span,
            input_data=params,
        )

    @classmethod
    def open(cls, path: str = None, primer_cost: PrimerCostModel = None):
        params = open_params(path)
        if not primer_cost:
            primer_cost = PrimerCostModel.open(path)
        return cls.from_json(params, primer_cost)

    def compute(self):
        gene_sizes = (
            self.syn[:, 0].reshape(-1, 1)[:: self.step_size, :].astype(np.int32)
        )
        gene_costs = self.syn[:, 1][gene_sizes]
        gene_times = self.syn[:, 2][gene_sizes]

        left_span = np.arange(self.lspanrange[0], self.lspanrange[1], self.step_size)[
            ..., np.newaxis, np.newaxis
        ]

        for i, j in [(0, 0), (1, 0), (0, 1), (1, 1)]:
            jxn = self._compute(gene_costs, gene_sizes, gene_times, i, j, left_span)
            self.cost_dict[(i, j)] = jxn

    def _compute(
        self,
        gene_costs,
        gene_sizes,
        gene_times,
        i: Union[bool, int],
        j: Union[bool, int],
        left_span,
    ):
        # extension conditions, idk
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
        _gtimes = syn_time_cost[idx[1]]
        gene_df = NumpyDataFrame(
            dict(
                cost=_gcosts,
                material=_gcosts,
                time=_gtimes,
                efficiency=np.ones(idx[0].shape[0]),
                size=gene_sizes[idx[1]],
            ),
            apply=np.squeeze,
        )

        flat_left_jxn = left_jxn[idx[2]].apply(np.squeeze)
        flat_right_jxn = right_jxn[idx[2], idx[1], idx[0]]

        time = np.vstack(
            (
                flat_left_jxn.data["time"],
                flat_right_jxn.data["time"],
                gene_df.data["time"],
            )
        ).max(axis=0)

        gap_df = NumpyDataFrame(
            dict(
                span=_span,
                cost=syn_total_cost[idx],
                efficiency=syn_eff[idx],
                time=time,
                material=syn_material_cost[idx],
                lshift=left_span[idx[2]],
            ),
            apply=np.squeeze,
        )

        gap_df.update(flat_left_jxn.prefix("lprimer_"))
        gap_df.update(flat_right_jxn.prefix("rprimer_"))
        gap_df.update(gene_df.prefix("gene_"))

        return gap_df

    def cost(
        self, bp: Union[int, np.ndarray], ext: Tuple[int, int], invalidate=True
    ) -> NumpyDataFrame:
        return super().cost(bp, ext, invalidate=invalidate)

    def __call__(self, bp, ext):
        return self.cost(bp, ext)


class SpanCost(CostBuilder):
    def __init__(
        self, syn_cost, input_data: dict = None, override_span: np.ndarray = None
    ):
        self.syn_cost = syn_cost
        self.primer_cost = self.syn_cost.primer_cost
        x = [
            syn_cost.span.min(),
            syn_cost.span.max(),
            self.primer_cost.span.min(),
            self.primer_cost.span.max(),
        ]
        if override_span is not None:
            span = override_span
        else:
            span = np.arange(min(x), max(x))
        self.input_data = input_data
        super().__init__(span)
        self.compute()

    @classmethod
    def from_json(cls, params, override_span: np.ndarray = None):
        super().from_json(params)
        primer_cost_model = PrimerCostModel.from_json(params)
        synthesis_cost_model = SynthesisCostModel.from_json(
            params, primer_cost_model, override_span=override_span
        )
        return cls(
            syn_cost=synthesis_cost_model,
            input_data=params,
            override_span=override_span,
        )

    def compute(self):
        def choose(a, i):
            return np.choose(i, a)

        for ext in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            # numpy data frames for primer cost and syn cost over span
            df1 = self.primer_cost(self.span, ext)
            df2 = self.syn_cost(self.span, ext)

            # determine the indices of the min cost (0=primer, 1=syn)
            c1 = df1.data["cost"]
            c2 = df2.data["cost"]
            c3 = np.stack((c1, c2), axis=1)
            y = c3.argmin(axis=1)

            # select between primer_cost and syn_cost based on the min cost
            df4 = NumpyDataFrame.group_apply(
                (df1, df2), choose, i=y, _fill_value=np.nan
            )
            self.cost_dict[ext] = df4

    def cost(
        self, bp: Union[int, np.ndarray], ext: Tuple[int, int], invalidate=True
    ) -> NumpyDataFrame:
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
