import hashlib
import os
from glob import glob

import pytest

from dasi import cost
from dasi.log import logger

logger.set_level("DEBUG")


def hashfiles(files, hash="sha1", hashfunc=None):

    if not hashfunc:

        def hashfunc(string):
            return getattr(hashlib, hash)(string.encode("utf-8")).hexdigest()

    contents = ""
    sorted_path = sorted(files)
    for file in sorted_path:
        with open(file, "r") as f:
            contents += hashfunc(f.read())
    return hashfunc(contents)


def cost_module_checksum():
    cost_dir = os.path.dirname(cost.__file__)
    cost_files = sorted(
        glob(os.path.join(cost_dir, "*.py")) + glob(os.path.join(cost_dir, "*.json"))
    )
    return hashfiles(cost_files)


here = os.path.abspath(os.path.dirname(__file__))
path = os.path.join(here, "span_cost.b")
checksum_path = os.path.join(here, "cost.checksum")


def cached(path, save_func, load_func, checksum_path):
    different_checksum = True
    checksum = cost_module_checksum()
    if os.path.isfile(checksum_path):
        with open(checksum_path, "r") as f:
            stored_checksum = f.read().strip()
            if stored_checksum == checksum:
                different_checksum = False
            logger.debug("Stored checksum: {}".format(stored_checksum))
    logger.debug("Checksum: {}".format(checksum))
    if different_checksum or not os.path.isfile(path):
        logger.debug("Using default params")
        model = save_func(path)
        # span_cost = cost.SpanCost.open()
        # span_cost.dump(path)
        with open(checksum_path, "w") as f:
            f.write(checksum)
        stat = os.stat(path)
        logger.debug("Wrote {} bytes".format(stat.st_size))
    else:
        logger.debug("Loading {}".format(path))
        stat = os.stat(path)
        logger.debug("Loaded {} bytes".format(stat.st_size))
        model = load_func(path)
        # span_cost = cost.SpanCost.load(path)
    return model


@pytest.fixture(scope="function")
def cached_span_cost():
    def load_span_cost(path):
        span_cost = cost.SpanCost.load(path)
        return span_cost

    def save_span_cost(path):
        span_cost = cost.SpanCost.open()
        span_cost.dump(path)
        return span_cost

    return cached(load_func=load_span_cost, save_func=load_span_cost)


@pytest.fixture(scope="module")
def primer_cost():
    return cost.PrimerCostModel.open()


@pytest.fixture(scope="module")
def syn_cost(primer_cost):
    return cost.SynthesisCostModel.open(primer_cost=primer_cost)


@pytest.fixture(scope="module")
def span_cost(primer_cost, syn_cost):
    return cost.SpanCost(syn_cost)


import dill


def test_dump_load_primer_cost(primer_cost):
    dill.loads(dill.dumps(primer_cost))


def test_dump_load_syn_cost(syn_cost):
    dill.loads(dill.dumps(syn_cost))
