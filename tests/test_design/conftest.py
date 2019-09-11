import pytest
from dasi.cost import SpanCost
import os
from dasi.log import logger

here = os.path.abspath(os.path.dirname(__file__))
do_save = True


@pytest.fixture(scope='module')
def span_cost():
    """Saves the span cost as bytes; reloads when called."""
    path = os.path.join(here, 'span_cost.b')
    if do_save and os.path.isfile(path):
        with logger.timeit("INFO", "loading bytes"):
            print("Loading file: {}".format(path))
            span_cost = SpanCost.load(path)
    else:
        span_cost = SpanCost.default()
        if do_save:
            with logger.timeit("INFO", "saving bytes"):
                print("Saving file: {}".format(path))
                span_cost.dump(path)
    return span_cost