import pytest
import dill
from dasi.cost import SpanCost
import os

here = os.path.abspath(os.path.dirname(__file__))
do_save = True


@pytest.fixture(scope='module')
def span_cost():
    path = os.path.join(here, 'span_cost.pkl')
    if do_save and os.path.isfile(path):
        with open(path, 'rb') as f:
            print("Loading file: {}".format(path))
            span_cost = dill.load(f)
    else:
        span_cost = SpanCost.default()
        if do_save:
            with open(path, 'wb') as f:
                dill.dump(span_cost, f)
    return span_cost