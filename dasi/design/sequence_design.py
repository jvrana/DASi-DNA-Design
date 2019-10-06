"""sequence_design.py.

A module for converting DASi designs into actual DNA sequences.
"""
from typing import Any
from typing import Dict
from typing import Tuple
from typing import Union

import networkx as nx
import primer3plus

from dasi.cost import SpanCost
from dasi.design import Assembly
from dasi.design import AssemblyNode
from dasi.utils import NumpyDataFrame
from dasi.utils import Region


def design_primers(
    template: str, region: Region, lseq: Union[None, str], rseq: Union[None, str]
) -> Tuple[Dict[int, dict], Dict[str, Any]]:
    """Design primers flanking the specified.

    :class:`Region.<dasi.utils.Region>`. If the region is cyclic and spans the
    origin, this method will handle the appropriate manipulations to design
    primers around the origin and restore the locations of the resulting primer
    pairs.

    :param template: the template string to design primers
    :param region: region specified to design primers around. Regions are exclusive at
                    their end points (`.b` parameter)
    :param lseq: optionally provided left sequence
    :param rseq: optionally provided right sequence
    :return: tuple of pairs and the 'explain' dictionary.
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
    design.presets.product_size((len(region), len(region)))
    design.PRIMER_PICK_ANYWAY = True
    design.presets.use_overhangs()
    design.presets.long_ok()
    print(
        primer3plus.utils.anneal(
            design.SEQUENCE_TEMPLATE.value, [design.SEQUENCE_PRIMER_REVCOMP.value]
        )
    )
    print(design.SEQUENCE_PRIMER.value)
    print(design.SEQUENCE_PRIMER_REVCOMP.value)
    print(design.SEQUENCE_INCLUDED_REGION.value)
    print(design.SEQUENCE_TEMPLATE.value)
    pairs, explain = design.run_and_optimize(15)
    if index is not None:
        for pair in pairs.values():
            loc = pair["LEFT"]["location"]
            pair["LEFT"]["location"] = (index[loc[0]], loc[1])

            loc = pair["RIGHT"]["location"]
            pair["RIGHT"]["location"] = (index[loc[0]], loc[1])
    return pairs, explain


def edata_to_npdf(edata: dict, span_cost: SpanCost) -> NumpyDataFrame:
    return span_cost.cost(edata["span"], edata["type_def"])


def design_pcr_product_primers(
    assembly: Assembly, n1: AssemblyNode, n2: AssemblyNode, span_cost: SpanCost
):
    pass
    # graph = assembly.graph
    # successors = graph.successors(n1)
    # predecessors = graph.predecessors(n2)
    #
    # edge = graph[n1][n2]
    # sedge = graph[n2][successors[0]]
    # pedge = graph[predecessors[0]][n1]
    #
    # # the left and right extensions are determined by the pred and succ edges.
    # left = edata_to_npdf(pedge, span_cost)
    # right = edata_to_npdf(sedge, span_cost)
    # print(edge)
    # print(left)
    # print(right)
