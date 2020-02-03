r"""
Schemas (:mod:`dasi.schemas`)
=============================

Validated JSON Schemas (https://json-schema.org/)

Check out the :ref:`JSON Schema page <json_schema>` for more information.
"""
import json
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Dict

import jsonschema

schema_directory = dirname(abspath(__file__))  #: path to directories containing schemas


def _load(filename: str) -> Dict:
    with open(join(schema_directory, filename)) as f:
        return json.load(f)


class Schemas:
    """Class that contains JSON schemas."""

    output_schema = _load("output_schema.json")
    cost_parameters_schema = _load("parameter_json_schema.json")


def validate_with_schema(
    data: Dict, schema: Dict, do_raise: bool = True, reraise_as: Exception = None
):
    """Validate JSON data using a JSON Schema."""
    try:
        jsonschema.validate(data, schema)
    except jsonschema.ValidationError as e:
        if do_raise:
            if reraise_as:
                raise reraise_as(str(e)) from e
            else:
                raise e
        else:
            return False
    return True
