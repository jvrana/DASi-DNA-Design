:tocdepth: 5

.. _api:

*************
API Reference
*************

.. module:: dasi

.. _api_design:

Design
======

.. automodule:: dasi.design
    :members:

.. autoclass:: dasi.design.Design
    :members:
    :inherited-members:

.. autoclass:: dasi.design.AssemblyGraphBuilder
    :members:
    :inherited-members:

.. autoclass:: dasi.design.Assembly
    :members:
    :inherited-members:

.. autoclass:: dasi.design.DesignResult
    :members:
    :inherited-members:

.. autoclass:: dasi.LibraryDesign
    :members:
    :inherited-members:

.. automodule:: dasi.design.design_algorithms
    :members:
    :private-members:

.. autoclass:: dasi.design.plotter.Plotter
    :members:

.. automodule:: dasi.design.sequence_design
    :members:


Alignments
==========


.. automodule:: dasi.alignments
    :members:
    :imported-members:


Cost Calculations
=================

.. autoclass:: dasi.cost.CostBuilder
    :members:
    :inherited-members:

.. autoclass:: dasi.cost.PrimerCostBuilder
    :members:

.. autoclass:: dasi.cost.SynthesisCostBuilder
    :members:

.. autoclass:: dasi.cost.SpanCost
    :members:

Parameters
----------

.. autoclass:: dasi.cost.Globals
    :members:

.. autoclass:: dasi.cost.PrimerParams
    :members:

.. autoclass:: dasi.cost.SynthesisParams
    :members:

Utilities
=========

.. automodule:: dasi.utils
    :members:

Cost utilities
--------------

.. automodule:: dasi.cost.utils
    :members:
    :imported-members:

Networkx utilities
------------------

.. automodule:: dasi.utils.networkx
    :members:
    :imported-members:
