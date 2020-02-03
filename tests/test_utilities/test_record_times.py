from dasi.utils import log_times


class Foo:
    def __init__(self):
        self._method_run_times = {}

    @log_times()
    def bar(self):
        return 1


def test_log_times():
    foo = Foo()
    foo.bar()
    print(foo._method_run_times)
