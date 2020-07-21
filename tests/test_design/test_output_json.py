import json

import pytest

from dasi import Design
from dasi import LibraryDesign
from dasi.design.output import dasi_design_to_output_json
from dasi.design.output import validate_output


@pytest.fixture(
    scope="module", params=[Design, LibraryDesign], ids=["Design", "LibraryDesign"]
)
def design(request):
    _Design = request.param
    design = _Design.fake(
        n_designs=3, n_linear_seqs=500, n_primers_from_templates=50, shared_length=0
    )
    design.run()
    return design


class TestOutput:
    def test_output(self, design):
        out = dasi_design_to_output_json(design)
        print(json.dumps(out, indent=2))

    def test_validate_output(self, design):
        out = dasi_design_to_output_json(design)
        validate_output(out)

    def test_to_out_json(self, design):
        out = design.out(elim_extra_reactions=True)
        assert out

    def test_to_out_json_dumps(self, design):
        out = design.out()
        assert json.dumps(out)
