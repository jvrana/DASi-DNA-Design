from bisect import bisect_left, bisect_right
from collections import Sized
from more_itertools import unique_everseen, flatten
from itertools import islice


def argsort(a):
    x = list(enumerate(a))
    x.sort(key=lambda b: b[1])
    return [_x[0] for _x in x]


def argsortandval(a):
    x = list(enumerate(a))
    x.sort(key=lambda b: b[1])
    x = tuple(zip(*[_x for _x in x]))
    if x:
        return x
    return [], []


class SpanDictionary(Sized):

    __slots__ = ["_values", "_mins", "_maxs"]

    def __init__(self, data=None):
        self._values = []
        self._mins = []
        self._maxs = []
        if data:
            for k, v in data.items():
                self[k[0], k[1]] = v

    def keys(self):
        return self.ranges()

    def values(self):
        return self._values[:]

    def ranges(self, data=False):
        if not data:
            return zip(self._mins, self._maxs)
        return zip(self._mins, self._maxs, self.values())

    def items(self, step=1):
        if step != 1:
            a = islice(self, None, None, step)
        else:
            a = self
        for k in a:
            yield k, self[k]

    def __setitem__(self, key, val):
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("Must set using [mn, mx]")
        if key[0] > key[1]:
            raise ValueError("Min must be creater than max when setting indices.")
        i = bisect_left(self._mins, key[0])
        self._values.insert(i, val)
        self._mins.insert(i, key[0])
        self._maxs.insert(i, key[1])

    def _get_indices(self, key):
        i = bisect_right(self._mins, key)

        indices = []
        for j, mx in enumerate(self._maxs[:i]):
            if key < mx:
                indices.append(j)
        return indices

    def __getitem__(self, key):
        i = bisect_right(self._mins, key)
        values = []
        for j, mx in enumerate(self._maxs[:i]):
            if key < mx:
                values.append(self._values[j])
        return values

    def to_dict(self):
        return dict(zip(self.ranges(), self._values))

    def min(self):
        return min(self._mins)

    def max(self):
        return max(self._maxs)

    def __contains__(self, item):
        if self._get_indices(item):
            return True
        return False

    def __iter__(self):
        return unique_everseen(flatten(range(a[0], a[1]) for a in self.ranges()))

    def __len__(self):
        return len(self._values)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{}({})".format(self.__class__.__name__, str(self.to_dict()))
