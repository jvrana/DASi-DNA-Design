r"""
.. _cost_model:

Cost Model (:mod:`dasi.cost`)
=============================

.. currentmodule:: dasi.cost

This module provide cost calculations for molecular assemblies.

Default cost parameters can be opened as follows:

.. code-block::

    cost_model = SpanCost.open()

    # optionally use custom parameters
    cost_model = SpanCost.open("my_parameters.json")

Take a look at the :ref:`parameters schema <cost_schema>` for how to build the
parameter json.

Utilities
---------

.. autosummary::
    :toctree: generated/

    utils
"""
from dasi.cost.cached_span_cost import cached_span_cost
from dasi.cost.span_cost import open_params
from dasi.cost.span_cost import PrimerCostModel
from dasi.cost.span_cost import SpanCost
from dasi.cost.span_cost import SynthesisCostModel
from dasi.cost.span_cost import validate_params
