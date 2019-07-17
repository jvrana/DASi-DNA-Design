from collections.abc import Container, Iterable, Sized
from itertools import chain
from copy import copy


class SpanError(Exception):
    pass


# TODO: move region code to lobio


class Span(Container, Iterable, Sized):
    __slots__ = ["a", "b", "context_length", "cyclic", "index"]

    def __init__(self, a, b, l, cyclic=False, index=0, allow_wrap=False):
        if a > b and not cyclic:
            raise IndexError(
                "Start {} cannot exceed end {} for linear spans".format(a, b)
            )

        self.index = index
        self.context_length = l
        self.cyclic = cyclic
        if not cyclic or not allow_wrap:
            bounds = self.bounds()
            if not bounds[0] <= a <= bounds[1]:
                raise IndexError(
                    "Start {} must be within ({}, {})".format(a, index, index + l)
                )
            if not bounds[0] <= b <= bounds[1]:
                raise IndexError(
                    "End {} must be within ({}, {})".format(b, index, index + l)
                )
        if a - index > l or a - index < 0:
            self.a = self.t(a)
        else:
            self.a = a
        if b - index > l or b - index < 0:
            self.b = self.t(b)
        else:
            self.b = b

    def bounds(self):
        """Return the bounds (end exclusive)"""
        return self.index, self.context_length + self.index

    def t(self, p):
        """Translate an index to a valid index"""
        if p >= self.bounds()[1] or p < self.bounds()[0]:
            if not self.cyclic:
                raise IndexError(
                    "Position {} outside of linear bounds {}".format(p, self.bounds())
                )
        return (p % self.context_length) + self.bounds()[0]

    def ranges(self):
        if self.cyclic and self.a > self.b:
            return [(self.a, self.bounds()[1]), (self.bounds()[0], self.b)]
        else:
            return [(self.a, self.b)]

    def new(self, a, b, allow_wrap=False):
        return self.__class__(
            a,
            b,
            self.context_length,
            self.cyclic,
            index=self.index,
            allow_wrap=allow_wrap,
        )

    def sub(self, a, b):
        if b is not None and a > b and not self.cyclic:
            raise ValueError(
                "Start {} cannot exceed end {} for linear spans".format(a, b)
            )
        copied = copy(self)
        copied.a = a
        copied.b = b
        return copied

    def same_context(self, other):
        return (
            other.context_length == self.context_length and self.cyclic == other.cyclic
        )

    def force_context(self, other):
        if not self.same_context(other):
            raise SpanError("Cannot compare with different contexts")

    def overlaps_with(self, other):
        self.force_context(other)
        # if other in self:
        #     return True
        # elif self in other:
        #     return True
        if (
            other.a in self
            or other.b - 1 in self
            or self.a in other
            or self.b - 1 in other
        ):
            return True
        return False

    def differences(self, other):
        """Return a tuple of differences between this span and the other span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return (self.sub(self.a, other.a),)
        elif other.b in self and other.a not in self:
            return (self.sub(other.b, self.b),)
        if other in self:
            return self.sub(self.a, other.a), self.sub(other.b, self.b)
        elif self in other:
            return (self.sub(self.a, self.a),)
        else:
            return ((self[:]),)

    def intersection(self, other):
        """Return the span inersection between this span and the other span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return self.sub(other.a, self.b)
        elif other.b in self and other.a not in self:
            return self.sub(self.a, other.b)
        if other in self:
            return other[:]
        elif self in other:
            return self[:]

    def consecutive(self, other):
        self.force_context(other)
        try:
            return self.b == other.t(other.a)
        except IndexError:
            return False

    def connecting_span(self, other):
        """
        Return the span that connects the two spans. Returns None

        :param other:
        :return:
        """
        self.force_context(other)
        if self.cyclic and self.a == other.a and self.b == other.b:
            return self.invert()[0]
        if self.consecutive(other):
            return self[self.b, self.b]
        elif self.overlaps_with(other):
            return None
        else:
            if self.b > other.a and not self.cyclic:
                return None
            return self[self.b, other.a]

    # def union(self, other):
    #     self.force_context(other)
    #     if other in self:
    #         return self[:]
    #     elif self in other:
    #         return other[:]
    #     elif other.a in self and not other.t(other.b - 1) in self:
    #         if self.a == other.b:
    #             return self.new(self.a, None)
    #         return self.new(self.a, other.b)
    #     elif other.t(other.b - 1) in self and other.a not in self:
    #         if other.a == self.b:
    #             return self.new(self.a, None)
    #         return self.new(other.a, self.b)

    # def __ge__(self, other):
    #     self.force_context(other)
    #     return self.a >= other.a

    #     def __lt__(self, other):
    #         return self.a < other.a

    #     def __gt__(self, other):
    #         return self.a > other.a

    #     def __ge__(self, other):
    #         return self.a >= other.a

    # def __invert__(self):
    #     if self.a > self.b:
    #         if self.cyclic:
    #             return self[self.b+1, self.a-1],
    #         else:
    # return

    def invert(self):
        if self.cyclic:
            return (self[self.b, self.a],)
        else:
            return self[:, self.a], self[self.b, :]

    def __eq__(self, other):
        return self.same_context(other) and self.a == other.a and self.b == other.b

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __contains__(self, other):
        if isinstance(other, int):
            for r in self.ranges():
                if r[0] <= other < r[1]:
                    return True
            return False
        elif issubclass(type(other), Span):
            if not self.same_context(other):
                return False
            return (
                other.a in self
                and other.t(other.b - 1) in self
                and len(other) <= len(self)
            )

    def __len__(self):
        return sum([r[1] - r[0] for r in self.ranges()])

    def __iter__(self):
        for i in chain(*[range(*x) for x in self.ranges()]):
            yield i

    def __invert__(self):
        return self.invert()

    def __getitem__(self, val):
        if isinstance(val, int):
            if not self.cyclic and (val >= len(self) or val < -len(self)):
                raise IndexError(
                    "Index '{}' outside of linear span with length of {}".format(
                        val, len(self)
                    )
                )
            else:
                if val < 0:
                    return self.t(val + self.b) - self.index
                else:
                    return self.t(val + self.a) - self.index
        elif issubclass(type(val), slice):
            if val.step:
                raise ValueError(
                    "{} slicing does not support step.".format(self.__class__.__name__)
                )
            if val.start is None:
                i = self.a
            else:
                i = self[val.start]

            if val.stop is None:
                j = self.b
            else:
                j = self[val.stop]
            return self.sub(i, j)
        elif isinstance(val, tuple):
            if len(val) > 2:
                raise ValueError(
                    "{} -- copying does only supports (start, stop)".format(val)
                )
            val = list(val)
            for i in [0, 1]:
                if val[i] == slice(None, None, None):
                    val[i] = None
            if val == (None, None):
                val = self.bounds()
            elif val[0] is None:
                val = (self.bounds()[0], val[1])
            elif val[1] is None:
                val = (val[0], self.bounds()[1])
            return self.new(*val)
        else:
            raise ValueError("indexing does not support {}".format(type(val)))

    def __repr__(self):
        return "<{}={} {} {} start={}>".format(
            self.__class__.__name__,
            id(self),
            (self.a, self.b, self.context_length),
            self.cyclic,
            self.index,
        )

    def __str__(self):
        return "<{} {} {} start={}>".format(
            self.__class__.__name__,
            (self.a, self.b, self.context_length),
            self.cyclic,
            self.index,
        )


"""
Project: jdna
File: region
Author: Justin
Date: 2/21/17

Description: Basic functionality for defining regions of linear, circularized, or reversed
regions of a sequence.

"""


class Direction(object):
    FORWARD = 1
    REVERSE = -1
    BOTH = 0


class Region(Span):
    __slots__ = ["name", "id", "direction"]

    FORWARD = Direction.FORWARD
    REVERSE = Direction.REVERSE
    BOTH = Direction.BOTH

    def __init__(
        self,
        start,
        end,
        length,
        cyclic=False,
        index=0,
        direction=FORWARD,
        name=None,
        id=None,
        allow_wrap=False,
    ):
        self.name = name
        self.id = id
        assert direction in [self.FORWARD, self.REVERSE, self.BOTH]
        self.direction = direction
        super().__init__(
            start, end, length, cyclic=cyclic, index=index, allow_wrap=allow_wrap
        )

    @property
    def start(self):
        if self.direction == self.REVERSE:
            return self.b
        else:
            return self.a

    @property
    def end(self):
        if self.direction == self.REVERSE:
            return self.a
        else:
            return self.b
