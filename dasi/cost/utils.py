"""Cost utilities."""
import numpy as np


def square_broadcast(a, b, a_axis=0, b_axis=1):
    """Broadcast b to the a_axis of ndarray 'a' and perform a hstack."""
    c = np.broadcast_to(b, (a.shape[a_axis], b.shape[b_axis]))
    return np.hstack([a, c])


def df_to_np_ranged(min_col, max_col, df, cols=None, dtype=None):
    """Expand a pandas data frame, which has min and max span columns defined,
    into a numpy ndarray. Specific columns to include in the expansion can be
    defined in the `remaining_cols` argument, else all columns are used.

    For example:

    .. code-block::

      df = pd.DataFrame([
        [0, 10, 10], [10, 20, 100]
      ], columns=['min', 'max', 'x'])
      a = df_to_np_ranged('min', 'max', df, dtype=np.float64)
      assert a.shape == (20, 2)
      assert a[5] == [5., 10.]
      assert a[17] == [17., 100.]
    """
    vblocks = []
    mn, mx = min_col, max_col
    if cols is None:
        cols = [c for c in df.columns if c not in [mn, mx]]
    for _, row in df.iterrows():
        _a = np.arange(row[mn], row[mx], dtype=np.int32).reshape(-1, 1)
        cost_row = row[cols]
        vblocks.append(square_broadcast(_a, cost_row.to_numpy(dtype=dtype), 0, 0))

    a = np.vstack(vblocks)
    return a


def unshape(shape, axis):
    if isinstance(axis, int):
        axis = (axis,)
    return [s for i, s in enumerate(shape) if i not in axis]


def flatten_axis(a, axis=None, copy=True):
    """Flatten the array along the specified axes."""
    if axis is None:
        return a.flatten()
    if isinstance(axis, int):
        axis = (axis,)
    to_axis = tuple(range(len(axis)))
    g = np.prod(unshape(a.shape, axis))
    view = np.moveaxis(a, axis, to_axis).reshape(-1, g)
    if copy:
        return view.copy()
    else:
        return view


def duplicates(b, axis=None):
    """Find duplicates within the specified axes.

    If axis not provided, flatten array and return duplicates.
    """
    b = flatten_axis(b, axis, copy=False)
    s = np.sort(b)
    if axis:
        duplicate_list = []
        for x in b:
            duplicate_list.append(duplicates(x))
        return duplicate_list
    return np.unique(s[:-1][s[1:] == s[:-1]])


def lexargmin(a, axis):
    shape = a[0].shape
    reshaped = tuple([_a.reshape(shape[axis], -1) for _a in a])
    i = np.lexsort(reshaped, axis=1)
    i = i[:, :1]
    i = np.unravel_index(i, unshape(shape, axis))
    imin = (np.arange(shape[axis]), *[np.squeeze(_i) for _i in i])
    return imin


class Slicer:
    def __getitem__(self, item):
        return item


slicer = Slicer()
