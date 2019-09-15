from pyblast.utils import Span


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
        region_id=None,
        allow_wrap=True,
        does_wrap_origin=False
    ):
        self.name = name
        self.id = region_id
        assert direction in [self.FORWARD, self.REVERSE, self.BOTH]
        self.direction = direction
        super().__init__(
            start, end, length, cyclic=cyclic, index=index, allow_wrap=allow_wrap, does_wrap_origin=does_wrap_origin
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
