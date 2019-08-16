from .span import Span

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
        allow_wrap=True,
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
