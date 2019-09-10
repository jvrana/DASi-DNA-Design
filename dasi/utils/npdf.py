"""NumpyDataFrame"""
import pprint
from collections import OrderedDict
from collections.abc import Mapping, Iterable

import numpy as np
import pandas as pd


class Null(object):
    """Not None."""


class NumpyDataFrameException(Exception):
    pass


class NumpyDataFrameIndexer(Mapping):
    def __init__(self, df):
        self.df = df

    def __len__(self):
        return len(self.df.columns)

    def __iter__(self):
        for c in self.df.columns:
            yield c

    def __contains__(self, item):
        for c in self:
            if item == c:
                return True
        return False

    def __delitem__(self, key):
        del self.df.data[key]

    def __setitem__(self, col, val):
        self.df.data[col] = val
        self.df.validate()

    def __getitem__(self, cols):
        if isinstance(cols, str):
            cols = (cols,)
        elif not isinstance(cols, Iterable):
            cols = (cols,)
        data = {k: self.df.data[k] for k in self.df.data if k in cols}
        return self.df.__class__(data)


class NumpyDataFrame(Mapping):
    """The NumpyDataFrame is a class halfway between pandas and numpy. It has named columns, indexing, slicing,
    function applications, and mathematical operations. Unlike pandas however, it maintains the multi-dimensionality
    of underlying data (as np.ndarray), allowing broadcasting and complex indexing.

    Usage:

    **indexing and columns**

    All of the underlying arrays can be slices and indexed using the slicing operations.
    The native np.ndarray indexing is used for each column, meaning a *view* is returned
    with the same memory locations.

    .. code-block::

      df = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      df[0]
      df[0, 1]
      df[np.array([0, 1])]

    Columns can be selected and returned:

    .. code-block::

      df.col['A', 'B']    # return df with 'A' and 'B'
      df.col['A']         # return df with only 'A'
      print(list(df.col)) # return the column names
      print(df.columns)   # also returns the column names

    New columns can be added:

    .. code-block::

        df.col['A'] = np.arange(10)

    Columns can be deleted:

    .. code-block::

        del df.col['B']

    Add prefix or suffix to column names:

      df.prefix('prefix_', cols=['A']) # add prefix only to 'A', return new df
      df.suffix('__suffix')             # add suffix to all columns, return new df

    **apply**

    Functions can be apply to each column using `np.apply`.
    For example, the following applies `np.reshape` to all of the columns individually,
    returning a new df:

    .. code-block::

      df = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      df.apply(np.reshape, (-1, 1))

    Functions can be applied to tuples of all of the columns by `np.aggregate`.
    For example, the following stacks all of the columns horizontally, returning
    a new df:

    .. code-block::

      df = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      df.aggregate(np.hstack)

    Functions can be applied to grouped columns of multiple dataframes using `np.group_apply`.
    For example, the following applies stackes each column in each df horizontally:

    .. code-block::

      df1 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      df2 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      NumpyDataFrame.group_apply((df1, df2), np.hstack)

    In another example, we can apply np.divide to two dfs, using `expand=True` to expand
    the underlying arguments to properly run `np.divide`. The following two strategies are
    functionally equivalent:

    .. code-block::

      # strategy 1
      df1.group_apply((df1, df2), np.divide, expand=True)

      # strategy 2
      def div(a):
        return np.divide(a[0], a[1])
      df1.group_apply((df1, df2), div)

    **operations**

    Mathematical operations can be performed the same as you would np.ndarrays:

    .. code-block::

      	df1 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
        df2 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})

        df3 = df1 + df2  # add each element in each column
        df3 += 10 # add 10 to each element
        df3 *= 2. # multiply each element by 2.
        df3 * df3 # multiply each element in each column element wise
        df3 ** 2
        df3 ** df2

    **concatenations and appending**

    Dataframes can be concatenated together by the following:

    .. code-block::

      	df1 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
        df2 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
        NumpyDataFrame.concat((df1, df2))

    Dataframes with different column can be concatenated together by setting a fill value

    .. code-block::

        df1 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
        df2 = NumpyDataFrame({'A': np.arange(10), 'C': np.arange(10)})
        NumpyDataFrame.concat((df1, df2), fill_value=np.nan)

    **conversions**

    *to pandas*
    .. code-block::

        df1.to_df()

    *to numpy*
    .. code-block::

      df1.aggragate(np.hstack)

    *change dtype*

    .. code-block::
      df1.apply(np.astype, np.float64)

    """

    def __init__(self, data=None, apply=None):
        if data is None:
            data = {}
        self._data = OrderedDict()
        self._data.update({k: v for k, v in data.items() if v is not None})
        if apply:
            self.apply(apply, inplace=True)
        self.validate()

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, new):
        self._data = new
        self.validate()

    def validate(self):
        shapes = set(v.shape for v in self.data.values())
        if len(shapes) > 1:
            keys_and_shapes = {}
            for k, v in self.data.items():
                keys_and_shapes.setdefault(v.shape, list()).append(k)
            raise NumpyDataFrameException(
                "{} can only have one shape. Found the following shapes {}. If you want to sqeeze all"
                "of the data, set 'apply=np.squeeze'".format(
                    self.__class__, keys_and_shapes
                )
            )

    def prefix(self, s, cols=None, inplace=False):
        if cols is None:
            cols = self.columns
        return self.apply_to_col_names(lambda x: s + x, cols=cols, inplace=inplace)

    def suffix(self, s, cols=None, inplace=False):
        if cols is None:
            cols = self.columns
        return self.apply_to_col_names(lambda x: x + s, cols=cols, inplace=inplace)

    def apply_to_col_names(self, func, *args, cols=None, inplace=False, **kwargs):
        if cols is None:
            cols = self.columns
        data = {func(k, *args, **kwargs): v for k, v in self.data.items() if k in cols}
        if inplace:
            self.data = data
            return self
        return self.__class__(data)

    def aggregate(self, func, *args, cols=None, **kwargs):
        if cols is None:
            cols = self.columns
        collapsed = [self.data[c] for c in cols]
        return func(collapsed, *args, **kwargs)

    def apply(self, func, *args, astype=None, preprocess=None, inplace=False, **kwargs):
        data = {}
        for k, v in self.data.items():
            try:
                if preprocess:
                    data[k] = func(preprocess(v), *args, **kwargs)
                else:
                    data[k] = func(v, *args, **kwargs)
            except Exception as e:
                raise NumpyDataFrameException(
                    "Could not apply '{}' because '{} {}'".format(
                        func.__name__, type(e), e
                    )
                ) from e
        if inplace:
            if astype is not None and astype is not self.__class__:
                raise NumpyDataFrameException(
                    "Cannot convert from {} to {} while inplace=True".format(
                        self.__class__, astype
                    )
                )
            self.data = data
            return self
        if astype is None:
            astype = self.__class__
        return astype(data)

    @classmethod
    def merge(cls, others):
        df = cls()
        for a in others:
            df.update(a)
        return df

    @classmethod
    def concat(cls, others, fill_value=Null):
        return cls.group_apply(others, np.hstack, _fill_value=fill_value)

    def append(self, other):
        self.group_apply((other,), np.hstack)

    def fill_missing(self, cols, value):
        for c in cols:
            if c not in self.col:
                self.col[c] = np.array([value for _ in range(len(self))])

    @classmethod
    def group_apply(cls, others, func, *args, expand=False, _fill_value=Null, **kwargs):
        d = {}
        other_cols = set(tuple(sorted(o.columns)) for o in others)
        if len(other_cols) > 1:
            if _fill_value is Null:
                raise NumpyDataFrameException(
                    "Cannot apply to group. Different columns found: {}".format(
                        other_cols
                    )
                )
            else:
                all_cols = []
                for o in others:
                    all_cols += list(o.columns)
                all_cols = sorted(list(set(all_cols)))
                for o in others:
                    o.fill_missing(all_cols, _fill_value)
        for o in others:
            for k, v in o.data.items():
                d.setdefault(k, list()).append(v)
        if expand:
            data = {k: func(*v, *args, **kwargs) for k, v in d.items()}
        else:
            data = {k: func(v, *args, **kwargs) for k, v in d.items()}
        return cls(data)

    @classmethod
    def stack(cls, others, axis):
        return cls.group_apply(others, np.stack, axis=axis)

    @classmethod
    def hstack(cls, others):
        return cls.group_apply(others, np.hstack)

    @classmethod
    def vstack(cls, others):
        return cls.group_apply(others, np.vstack)

    @property
    def shape(self):
        return list(self.data.values())[0].shape

    def reshape(self, shape):
        return self.apply(np.reshape, shape)

    @property
    def columns(self):
        return tuple(self.data)

    @property
    def col(self):
        return NumpyDataFrameIndexer(self)

    def to_df(self, squeeze=True):
        if squeeze:
            return pd.DataFrame(self.apply(np.squeeze).data)
        return pd.DataFrame(self.data)

    def update(self, data, apply=None):
        if issubclass(type(data), NumpyDataFrame):
            data = data.data
        if issubclass(type(data), dict):
            self.data.update(data)
            if apply:
                self.apply(apply)
            self.validate()

    def items(self):
        return self.data.items()

    def copy(self):
        return self.apply(np.copy)

    def __getitem__(self, key):
        new = self.__class__(self.data)
        new.data = {k: v[key] for k, v in new.data.items()}
        return new

    def __setitem__(self, key, val):
        for v in self.data.values():
            v[key] = val

    def __len__(self):
        return len(list(self.data.values())[0])

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __add__(self, other):
        if not issubclass(type(other), NumpyDataFrame):
            return self.apply(np.sum, preprocess=lambda x: (x, other))
        return self.stack([self, other], axis=1).apply(np.sum, axis=1)

    def __mul__(self, other):
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.multiply, other)
        return self.group_apply([self, other], np.multiply, expand=True)

    def __truediv__(self, other):
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.divide, other)
        return self.group_apply([self, other], np.divide, expand=True)

    def __pow__(self, other):
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.power, other)
        return self.group_apply([self, other], np.power, expand=True)

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + -other

    def __str__(self):
        stacked = self.aggregate(np.stack, axis=1)
        return "{}\ncols={}\n{}".format(
            self.__class__, self.columns, pprint.pformat(stacked)
        )

    def __repr__(self):
        return str(self)
