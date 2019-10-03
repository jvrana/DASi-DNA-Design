"""sequence_design.py.

A module for converting DASi designs into actual DNA sequences.
"""
import operator
from itertools import accumulate

import primer3plus


def slice_with_region(iterable, region, func=operator.add):
    x = [iterable[s] for s in region.slices()]
    return accumulate(x, func)


def design_primers(template, region, lseq, rseq):
    design = primer3plus.new()
    design.presets.as_cloning_task()
    design.presets.template(template)
    design.presets.included((region.a, len(region)))
    # design.presets.product_size((len(region), len(region)))
    # design.presets.included([(region.a, len(region))])
    if lseq:
        design.presets.left_sequence(lseq)
    if rseq:
        design.presets.right_sequence(rseq)
    pairs, explain = design.run_and_optimize(5)
    return pairs, explain
