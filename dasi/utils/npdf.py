"""NumpyDataFrame."""
import pprint
from collections import OrderedDict
from collections.abc import Mapping as MappingABC
from typing import Any
from typing import Callable
from typing import Dict
from typing import Generator
from typing import ItemsView
from typing import Iterable
from typing import Iterable as IterableABC
from typing import List
from typing import Tuple
from typing import Union

import msgpack
import msgpack_numpy as m
import numpy as np
import pandas as pd

m.patch()


class Null:
    """Not None."""


class NumpyDataFrameException(Exception):
    """Generic exceptions for NumpyDataFrame."""


class NumpyDataFrame(MappingABC):
    """The NumpyDataFrame is a class halfway between pandas and numpy. It has
    named columns, indexing, slicing, function applications, and mathematical
    operations. Unlike pandas however, it maintains the multi-dimensionality of
    underlying data (as np.ndarray), allowing broadcasting and complex
    indexing.

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

    Functions can be applied to grouped columns of multiple dataframes using
    `np.group_apply`.
    For example, the following applies stackes each column in each df horizontally:

    .. code-block::

      df1 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      df2 = NumpyDataFrame({'A': np.arange(10), 'B': np.arange(10)})
      NumpyDataFrame.group_apply((df1, df2), np.hstack)

    In another example, we can apply np.divide to two dfs, using `expand=True` to expand
    the underlying arguments to properly run `np.divide`. The following two strategies
    are
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

    Dataframes with different column can be concatenated together by setting a fill
    value

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
        """Initializes a numpy data frame from a dict of string to np.ndarrays.
        The dict keys are representative of *column names* and the values are
        rows of that column. The shapes of each np.ndarray must be the same,
        else :class:`NumpyDataFrameException` is raised.

        :param data: A dict of string to np.ndarrays
        :param apply: Function to apply across the numpy data frame
        """
        if data is None:
            data = {}
        self._data = OrderedDict()
        self._data.update({k: v for k, v in data.items() if v is not None})
        if apply:
            self.apply(apply, inplace=True)
        self.validate()

    @property
    def data(self) -> Dict[str, np.ndarray]:
        """The underlying data dict of the dataframe."""
        return self._data

    @data.setter
    def data(self, new):
        """Set and validate the underlying data dict for the dataframe."""
        self._data = new
        self.validate()

    def validate(self):
        """Validate that the shapes of all of the np.ndarrays are the same."""
        shapes = {v.shape for v in self.data.values()}
        if len(shapes) > 1:
            keys_and_shapes = {}
            for k, v in self.data.items():
                keys_and_shapes.setdefault(v.shape, list()).append(k)
            raise NumpyDataFrameException(
                "{} can only have one shape. Found the following shapes {}. If you want"
                " to sqeeze all of the data, set 'apply=np.squeeze'".format(
                    self.__class__, keys_and_shapes
                )
            )
        for k, v in self.data.items():
            if not issubclass(type(v), np.ndarray):
                raise NumpyDataFrameException(
                    "{cls} only supports {typ}, not {x}".format(
                        cls=self.__class__.__name__, typ=np.ndarray, x=type(v)
                    )
                )

    def prefix(self, s: str, cols=None, inplace=False) -> "NumpyDataFrame":
        """Adds a prefix to all of the column names and returns a new
        dataframe."""
        if cols is None:
            cols = self.columns
        return self.apply_to_col_names(lambda x: s + x, cols=cols, inplace=inplace)

    def suffix(self, s: str, cols=None, inplace=False) -> "NumpyDataFrame":
        """Adds a prefix to all of the column names and returns a new
        dataframe."""
        if cols is None:
            cols = self.columns
        return self.apply_to_col_names(lambda x: x + s, cols=cols, inplace=inplace)

    def apply_to_col_names(
        self, func, *args, cols=None, inplace=False, **kwargs
    ) -> "NumpyDataFrame":
        """Apply a function to the column names and returns a new dataframe.

        :param func: the function to apply
        :param args: the additional function arguments
        :param cols: the columns to apply the function to. If False, all columns are
        used.
        :param inplace: if True, will apply the function to the current df and return
        the current df.
        :param kwargs: the additional function keyword arguments
        :return:
        """
        if cols is None:
            cols = self.columns
        data = {func(k, *args, **kwargs): v for k, v in self.data.items() if k in cols}
        if inplace:
            self.data = data
            return self
        return self.__class__(data)

    def aggregate(self, func, *args, cols=None, **kwargs) -> Any:
        """Group all of the np.ndarrays across all columns as a list and apply
        a function.

        :param func: the function to apply (e.g. np.hstack)
        :param args: the additional arguments of the function
        :param cols: the cols to apply the function to. If False, all columns are used.
        :param kwargs: the keyword arguments to apply to the function
        :return: the result of the function
        """
        if cols is None:
            cols = self.columns
        collapsed = [self.data[c] for c in cols]
        return func(collapsed, *args, **kwargs)

    def apply(
        self,
        func: Callable,
        *args,
        astype=None,
        preprocess=None,
        inplace=False,
        **kwargs,
    ) -> Any:
        """Apply a function to each np.ndarray.

        :param func: The function to appluy
        :param args: the additional arguments of the function
        :param astype: the type of data frame to return
        :param preprocess: preprocess function to apply to each np.ndarray before
        applying 'func'
        :param inplace: If True, will apply the function to the current df and return
        the current df.
        :param kwargs: the keyword arguments to apply to the function
        :return: a new dataframe
        """
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
    def merge(cls, others: Iterable["NumpyDataFrame"]) -> "NumpyDataFrame":
        """Merge many dfs into a single df."""
        df = cls()
        for a in others:
            df.update(a)
        return df

    @classmethod
    def concat(
        cls, others: Iterable["NumpyDataFrame"], fill_value=Null
    ) -> "NumpyDataFrame":
        """Concatenate several dfs into a single df."""
        return cls.group_apply(others, np.hstack, _fill_value=fill_value)

    def append(self, other: "NumpyDataFrame") -> "NumpyDataFrame":
        """Append the contents of the other df to this df."""
        self.group_apply((other,), np.hstack)
        return self

    def fill_value(self, cols: Iterable[str], value: Any) -> None:
        """Create new columns, if they are missing, and fill them with the
        specified value."""
        for c in cols:
            if c not in self.col:
                self.col[c] = np.array([value for _ in range(len(self))])

    @classmethod
    def group_apply(
        cls,
        others: Iterable["NumpyDataFrame"],
        func,
        *args,
        expand=False,
        _fill_value=Null,
        **kwargs,
    ) -> "NumpyDataFrame":
        """Groups np.arrays according to their column name for several
        dataframes (as a list) and applies a function to each group. Returns a
        new df with the results.

        :param others: iterable of dfs
        :param func: the function to apply
        :param args: additional arguments for the function
        :param expand: If true, the list of np.arrays will be expanded, as in
        `func(*list_of_arrs, ...)`
        :param _fill_value:
        :param kwargs: additional keyword arguments for the function
        :return: a new df
        """
        d = {}
        other_cols = {tuple(sorted(o.columns)) for o in others}
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
                    o.fill_value(all_cols, _fill_value)
        for o in others:
            for k, v in o.data.items():
                d.setdefault(k, list()).append(v)
        if expand:
            data = {k: func(*v, *args, **kwargs) for k, v in d.items()}
        else:
            data = {k: func(v, *args, **kwargs) for k, v in d.items()}
        return cls(data)

    @classmethod
    def stack(cls, others: Iterable["NumpyDataFrame"], axis):
        """Apply np.stack to each column for several dfs."""
        return cls.group_apply(others, np.stack, axis=axis)

    @classmethod
    def hstack(cls, others: Iterable["NumpyDataFrame"]):
        """Apply np.hstack to each column for several dfs."""
        return cls.group_apply(others, np.hstack)

    @classmethod
    def vstack(cls, others: List["NumpyDataFrame"]):
        """Apply np.vstack to each column for several dfs."""
        return cls.group_apply(others, np.vstack)

    @property
    def shape(self) -> Tuple[int, ...]:
        """Return the expected shape for the underlying np.ndarray.

        This is the shape of the array for any given column.
        """
        return list(self.data.values())[0].shape

    def reshape(self, shape) -> "NumpyDataFrame":
        """Reshape all arrays in the df."""
        return self.apply(np.reshape, shape)

    @property
    def columns(self) -> Tuple[str, ...]:
        """Return the column names."""
        return tuple(self.data)

    @property
    def col(self) -> "NumpyDataFrameIndexer":
        """Return the column indexer."""
        return NumpyDataFrameIndexer(self)

    def to_df(self, force=True) -> pd.DataFrame:
        """Force the NumpyDataFrame into a pandas.DataFrame."""
        print(self.shape)
        x = self
        if force:
            if len(x.shape) >= 2:
                x = x.apply(np.squeeze)
            elif len(x.shape) == 0:
                x = x.apply(np.expand_dims, axis=0)
        if len(x.shape) > 1:
            raise NumpyDataFrameException(
                "Unable to force dataframe of "
                "shape {} into a pandas.DataFrame".format(self.shape)
            )
        return pd.DataFrame(x.data)

    def update(self, data: Union["NumpyDataFrame", Dict[str, np.ndarray]], apply=None):
        """Update the df from a dict or another df."""
        if issubclass(type(data), NumpyDataFrame):
            data = data.data
        if issubclass(type(data), dict):
            self.data.update(data)
            if apply:
                self.apply(apply)
            self.validate()

    def items(self) -> ItemsView[str, np.ndarray]:
        """Iterate key: arr for the the underlying data dict."""
        return self.data.items()

    def copy(self) -> "NumpyDataFrame":
        """Copy the df."""
        return self.apply(np.copy)

    def dumps(self):
        """Use msgpack to dump df to a byte string."""
        return msgpack.dumps(self.data)

    def dump(self, f: str):
        """Dump byte repr of df to the specified path."""
        with open(f, "wb") as f:
            msgpack.dump(self.data, f)

    @classmethod
    def loads(cls, s: str) -> "NumpyDataFrame":
        """Use msgpack to load a df from a byte string."""
        data = msgpack.loads(s)
        data = {k.decode(): v for k, v in data.items()}
        return cls(data)

    @classmethod
    def load(cls, f: str) -> "NumpyDataFrame":
        """Load the byte repr of df from the specified path."""
        with open(f, "rb") as f:
            data = msgpack.load(f)
            data = {k.decode(): v for k, v in data.items()}
            return cls(data)

    def __getitem__(self, key: Union[int, slice, np.ndarray]) -> "NumpyDataFrame":
        new = self.__class__(self.data)
        new.data = {k: np.array(v[key]) for k, v in new.data.items()}
        return new

    def __setitem__(self, key: int, val: Any):
        for v in self.data.values():
            v[key] = val

    def __len__(self) -> int:
        return self.shape[0]

    def __iter__(self) -> Generator["NumpyDataFrame", None, None]:
        for i in range(len(self)):
            yield self[i]

    def __add__(self, other: "NumpyDataFrame") -> "NumpyDataFrame":
        if not issubclass(type(other), NumpyDataFrame):
            return self.apply(np.sum, preprocess=lambda x: (x, other))
        return self.stack([self, other], axis=1).apply(np.sum, axis=1)

    def __mul__(self, other: Union[int, float]) -> "NumpyDataFrame":
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.multiply, other)
        return self.group_apply([self, other], np.multiply, expand=True)

    def __truediv__(self, other: "NumpyDataFrame") -> "NumpyDataFrame":
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.divide, other)
        return self.group_apply([self, other], np.divide, expand=True)

    def __pow__(self, other: "NumpyDataFrame") -> "NumpyDataFrame":
        if not issubclass(type(other), NumpyDataFrame):
            other = np.array([other] * self.shape[0])
            return self.apply(np.power, other)
        return self.group_apply([self, other], np.power, expand=True)

    def __neg__(self) -> "NumpyDataFrame":
        return self * -1

    def __sub__(self, other: "NumpyDataFrame") -> "NumpyDataFrame":
        return self + -other

    def __str__(self) -> str:
        if self.shape == tuple():
            x = self.apply(np.expand_dims, axis=0)
            x = x.aggregate(np.stack, axis=0).flatten()
        else:
            x = self.aggregate(np.stack, axis=1)
        return "{}\ncols={}\n{}".format(self.__class__, self.columns, pprint.pformat(x))

    def __repr__(self) -> str:
        return "<{cls} shape={shape} cols={cols} dtype={dtype}>".format(
            cls=self.__class__.__name__,
            shape=self.shape,
            cols=self.columns,
            dtype=list(self.data.values())[0].dtype,
        )


class NumpyDataFrameIndexer(MappingABC):
    """The indexer for NumpyDataFrames."""

    def __init__(self, df: "NumpyDataFrame"):
        self.df = df

    def __len__(self) -> int:
        """Return number of columns."""
        return len(self.df.columns)

    def __iter__(self) -> Generator[str, None, None]:
        yield from self.df.columns

    def __contains__(self, item: str) -> bool:
        for c in self:
            if item == c:
                return True
        return False

    def __delitem__(self, key):
        del self.df.data[key]

    def __setitem__(self, col: str, val: np.ndarray):
        if not issubclass(type(val), np.ndarray):
            raise TypeError("Value must be a np.ndarray, not a {}".format(type(val)))
        self.df.data[col] = val
        self.df.validate()

    def __getitem__(self, cols: Union[str, Iterable[str]]) -> "NumpyDataFrame":
        if isinstance(cols, str):
            cols = (cols,)
        elif not isinstance(cols, IterableABC):
            cols = (cols,)
        data = {k: self.df.data[k] for k in self.df.data if k in cols}
        return self.df.__class__(data)
