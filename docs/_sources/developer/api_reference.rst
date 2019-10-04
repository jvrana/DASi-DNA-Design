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

.. automodule:: dasi.cost
    :members: Globals, dasi.cost.params.PrimerParams, SynthesisParams

Parameters
----------

.. automodule:: dasi.cost.params

.. autoclass:: dasi.cost.params.Globals
    :members:

.. autoclass:: dasi.cost.params.PrimerParams
    :members:

.. autoclass:: dasi.cost.params.SynthesisParams
    :members:

Utilities
=========

.. automodule:: dasi.utils
    :members:
    :imported-members:

Cost utilities
--------------

.. automodule:: dasi.cost.utils
    :members:
    :imported-members:
    :private-members:

Networkx utilities
------------------

.. automodule:: dasi.utils.networkx
    :members:
    :imported-members:
    :private-members:
