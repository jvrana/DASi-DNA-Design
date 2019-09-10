import pandas as pd
import numpy as np
from collections.abc import Mapping
import pprint


class Junction(Mapping):

    def __init__(self, span, cost, material, efficiency, **kwargs):
        self.data = {
            'span': span,
            'cost': cost,
            'material': material,
            'efficiency': efficiency
        }
        self.data.update(kwargs)
        self.data = {k: v for k, v in self.data.items()}

    @property
    def shape(self):
        return list(self.data.values())[0].shape

    def to_df(self):
        darr = {k: np.squeeze(v) for k, v in self.data.items() if not issubclass(type(v), Junction)}
        df = pd.DataFrame(darr)

        keys = []
        other_dfs = []

        for k, v in self.data.items():
            if issubclass(type(v), Junction):
                keys.append(k)
                other_dfs.append(v.to_df())

        if keys:
            concat_df = pd.concat([df] + other_dfs, keys=['default'] + keys, axis=1)
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


class PrimerJunction(Junction):

    def __init__(self, span, cost, material, efficiency, left_ext, right_ext):
        super().__init__(span, cost, material, efficiency, left_ext=left_ext, right_ext=right_ext)


class GeneJunction(Junction):
    pass


class GapJunction(Junction):

    def __init__(self, span, cost, material, efficiency, gene, lprimer, rprimer, lshift):
        super().__init__(span, cost, material, efficiency, gene=gene, lprimer=lprimer, rprimer=rprimer, lshift=lshift)