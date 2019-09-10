# @title Cost Containers
import pandas as pd
import numpy as np
from collections.abc import Mapping
import pprint


class NumpyDataFrame(Mapping):
    def __init__(self, data, apply=None):
        self.data = {}
        self.data.update({k: v for k, v in data.items() if v is not None})
        if apply:
            self.data = {k: apply(v) for k, v in self.data.items()}
        self.validate()

    @property
    def shape(self):
        return list(self.data.values())[0].shape

    def validate(self):
        shapes = set(v.shape for v in self.data.values())
        if len(shapes) > 1:
            keys_and_shapes = {}
            for k, v in self.data.items():
                keys_and_shapes.setdefault(v.shape, list()).append(k)
            raise Exception(
                "Data can only have one shape. Found the following shapes {}. If you want to sqeeze all"
                "of the data, set 'apply=np.squeeze'".format(keys_and_shapes)
            )

    @property
    def columns(self):
        return tuple(self.data)

    @classmethod
    def stack(cls, a, axis):
        # collect
        tuple_dict = {}
        for _a in a:
            for k, v in _a.data.items():
                tuple_dict.setdefault(k, list()).append(v)

        stacked_dict = {}
        for k, v in tuple_dict.items():
            if isinstance(v[0], np.ndarray):
                stacked_dict[k] = np.stack(v, axis=axis)
            elif isinstance(v[0], Junction):
                stacked_dict[k] = v[0].stack(v, axis=axis)
        return cls(**stacked_dict)

    def squeeze(self):
        data = {k: v.squeeze() for k, v in self.data.items()}
        return self.__class__(**data)

    def flatten(self):
        data = {k: v.flatten() for k, v in self.data.items()}
        return self.__class__(**data)

    def to_df(self):
        darr = {
            k: np.squeeze(v)
            for k, v in self.data.items()
            if not issubclass(type(v), Junction)
        }
        df = pd.DataFrame(darr)

        keys = []
        other_dfs = []

        for k, v in self.data.items():
            if issubclass(type(v), self.__class__):
                keys.append(k)
                other_dfs.append(v.to_df())

        if keys:
            concat_df = pd.concat([df] + other_dfs, keys=["default"] + keys, axis=1)
            return concat_df
        else:
            return df

    def __getitem__(self, key):
        new = self.__class__(**self.data)
        new.data = {k: v[key] for k, v in new.data.items()}
        return new

    def __len__(self):
        return len(list(self.data.values())[0])

    def __iter__(self):
        return None

    def items(self):
        return self.data.items()

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return pprint.pprint(self.data)


class Junction(NumpyDataFrame):
    def __init__(
        self, span, cost, material, efficiency, condition=None, apply=None, **kwargs
    ):
        data_dict = dict(
            span=span,
            cost=cost,
            material=material,
            efficiency=efficiency,
            condition=condition,
        )
        data_dict.update(**kwargs)
        super().__init__(data_dict, apply=apply)


class PrimerJunction(Junction):
    def __init__(
        self,
        span,
        cost,
        material,
        efficiency,
        left_ext,
        right_ext,
        condition=None,
        apply=None,
    ):
        super().__init__(
            span,
            cost,
            material,
            efficiency,
            left_ext=left_ext,
            right_ext=right_ext,
            condition=condition,
            apply=apply,
        )


class GeneJunction(Junction):
    def __init__(self, span, cost, material, efficiency, size, apply=None):
        super().__init__(span, cost, material, efficiency, size=size, apply=apply)


class GapJunction(Junction):
    def __init__(
        self,
        span,
        cost,
        material,
        efficiency,
        gene,
        lprimer,
        rprimer,
        lshift,
        condition=None,
        apply=None,
    ):
        super().__init__(
            span,
            cost,
            material,
            efficiency,
            gene=gene,
            lprimer=lprimer,
            rprimer=rprimer,
            lshift=lshift,
            condition=condition,
            apply=apply,
        )
