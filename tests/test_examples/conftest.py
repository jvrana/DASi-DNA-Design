import os
import pathlib

import pytest


def _pytest_only_subdir_items(config, items):
    rootdir = pathlib.Path(config.rootdir)
    here = os.path.abspath(os.path.dirname(__file__))
    conftest_dir = pathlib.Path(here).relative_to(rootdir)
    for item in items:
        rel_path = pathlib.Path(item.fspath).relative_to(rootdir)
        if conftest_dir in rel_path.parents:
            yield item


def pytest_collection_modifyitems(config, items):
    # python 3.4/3.5 compat: rootdir = pathlib.Path(str(config.rootdir))
    for item in _pytest_only_subdir_items(config, items):
        item.add_marker(pytest.mark.slowtest)
