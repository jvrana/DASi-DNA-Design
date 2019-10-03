"""sequence_design.py.

A module for converting DASi designs into actual DNA sequences.
"""
import operator
from itertools import accumulate
from typing import Any
from typing import Dict
from typing import Tuple
from typing import Union

import primer3plus

from dasi.utils import Region


def design_primers(
    template: str, region: Region, lseq: Union[None, str], rseq: Union[None, str]
) -> Tuple[Dict[int, dict], Dict[str, Any]]:
    """Design primers flanking the specified :class:`Region.

    <dasi.utils.Region>`. If the region is cyclic and spans the origin, this
    method will handle the appropriate manipulations to design primers around
    the origin and restore the locations of the resulting primer pairs.

    :param template: the template string to design primers
    :param region: region specified to design primers around. Regions are exclusive at
                    their end points (`.b` parameter)
    :param lseq: optional provided left sequence
    :param rseq: optional provided right sequence
    :return:
    """
    design = primer3plus.new()
    design.presets.as_cloning_task()

    if region.spans_origin():
        adjusted_template = region.get_slice(template) + region.invert()[0].get_slice(
            template
        )
        design.presets.template(adjusted_template)
        design.presets.included((0, len(region)))
        index = list(region) + list(region.invert()[0])
    else:
        design.presets.template(template)
        design.presets.included((region.a, len(region)))
        index = None
    if lseq:
        design.presets.left_sequence(lseq)
    if rseq:
        design.presets.right_sequence(rseq)
    pairs, explain = design.run_and_optimize(5)
    if index is not None:
        for pair in pairs.values():
            loc = pair["LEFT"]["location"]
            pair["LEFT"]["location"] = (index[loc[0]], loc[1])

            loc = pair["RIGHT"]["location"]
            pair["RIGHT"]["location"] = (index[loc[0]], loc[1])
    return pairs, explain
