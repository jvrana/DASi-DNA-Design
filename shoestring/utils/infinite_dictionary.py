from collections import OrderedDict
from bisect import bisect_left
from sortedcontainers import SortedDict


class InfiniteDict(SortedDict):
    """
    An infinite sorted dictionary. If a key is not found, the value 'below' it is
    returned. If 'ceil=True', the value 'above' it is returned. In cases outside of
    bounds, the values at the lowest or highest bounds are returned
    """

    def __init__(self, *args, ceil=False, **kwargs):
        self._ceil = ceil
        super().__init__(*args, **kwargs)

    def __getitem__(self, k):
        if k in self:
            return super().__getitem__(k)
        else:
            i = self.bisect_left(k)
            if k >= self.keys()[-1]:
                return self.values()[-1]
            if self._ceil:
                return self.values()[i]
            else:

                return self.values()[max(i - 1, 0)]
