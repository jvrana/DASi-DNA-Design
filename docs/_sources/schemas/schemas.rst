.. _json_schema:

JSON Schemas
============

JSON inputs and schemas for DASi.

.. _cost_schema:

Default Cost Parameters
-----------------------

.. literalinclude:: cost_default.json
   :linenos:
   :language: json

Cost Parameter Schema
---------------------

The JSON schema for the :ref:`Cost Model <cost_model>` input. Initialize
the cost model using:

.. code-block:: python

    from dasi.cost import SpanCost
    # open custom parameters
    cost_model = SpanCost.open('my example.json')

    # open default parameters
    cost_model = SpanCost.open()

.. jsonschema:: cost_schema.json

Design Output Schema
--------------------

.. literalinclude:: example_out.json
   :linenos:
   :language: json

.. jsonschema:: output_schema.json
