import hashlib
import os
from glob import glob
from os.path import abspath
from os.path import dirname
from os.path import isfile
from os.path import join

from dasi.cost import span_cost
from dasi.cost.span_cost import SpanCost
from dasi.log import logger

here = abspath(dirname(__file__))
DEFAULT_COST_FILEPATH = join(here, "span_cost.b")
DEFAULT_COST_CHECKSUM = join(here, "span_cost.checksum")


def hashfiles(files, hash="sha1", hashfunc=None):
    if not hashfunc:

        def hashfunc(string):
            return getattr(hashlib, hash)(string.encode("utf-8")).hexdigest()

    contents = ""
    sorted_path = sorted(files)
    for file in sorted_path:
        with open(file) as f:
            contents += hashfunc(f.read())
    return hashfunc(contents)


def cost_module_checksum():
    cost_dir = dirname(span_cost.__file__)
    cost_files = sorted(glob(join(cost_dir, "*.py")) + glob(join(cost_dir, "*.json")))
    return hashfiles(cost_files)


def cached(path, save_func, load_func, checksum_path, logger=None):
    different_checksum = True
    checksum = cost_module_checksum()
    if os.path.isfile(checksum_path):
        with open(checksum_path) as f:
            stored_checksum = f.read().strip()
            if stored_checksum == checksum:
                different_checksum = False
            if logger:
                logger.debug("Stored checksum: {}".format(stored_checksum))
    if logger:
        logger.debug("Checksum: {}".format(checksum))
    if different_checksum or not isfile(path):
        if logger:
            logger.debug("Using default params")
        model = save_func(path)
        with open(checksum_path, "w") as f:
            f.write(checksum)
        stat = os.stat(path)
        if logger:
            logger.debug("Wrote {} bytes".format(stat.st_size))
    else:
        if logger:
            logger.debug("Loading {}".format(path))
        stat = os.stat(path)
        if logger:
            logger.debug("Loaded {} bytes".format(stat.st_size))
        model = load_func(path)
    return model


def cached_span_cost(
    cost_filepath=DEFAULT_COST_FILEPATH, cost_checksum_filepath=DEFAULT_COST_CHECKSUM
):
    """This will check the checksum of the cost module against the last
    checksum. If checksums are the same, the span cost will be loaded. Else,
    span_cost will be created from default parameters and saved with the cost
    module's checksum.

    :param cost_filepath: path of the span_cost
    :param cost_checksum_filepath: path of the checksum
    :return: SpanCost
    """

    def load_span_cost(path):
        span_cost = SpanCost.load(path)
        return span_cost

    def save_span_cost(path):
        span_cost = SpanCost.open()
        span_cost.dump(path)
        return span_cost

    return cached(
        cost_filepath,
        load_func=load_span_cost,
        save_func=save_span_cost,
        checksum_path=cost_checksum_filepath,
        logger=logger,
    )
