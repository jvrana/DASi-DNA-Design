"""sequence_design.py.

Methods for designing DNA sequences from the result of the Design module
"""
from typing import Any
from typing import Dict
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
import primer3plus
from Bio.SeqRecord import SeqRecord
from primer3plus.utils import reverse_complement as rc

from dasi.alignments import AlignmentGroup
from dasi.cost import SpanCost
from dasi.design.assembly import Assembly
from dasi.design.assembly import AssemblyNode
from dasi.exceptions import DasiNoPrimerPairsException
from dasi.exceptions import DasiSequenceDesignException
from dasi.utils import NumpyDataFrame
from dasi.utils import Region


def design_primers(
    template: str,
    region: Region,
    lseq: Union[None, str],
    rseq: Union[None, str],
    left_overhang: Union[None, str] = None,
    right_overhang: Union[None, str] = None,
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
    :param loverhang: optionally provided left overhang sequence of the primer
    :param roverhang: optionally provided right overhang sequence of the primer
    :return: tuple of pairs and the 'explain' dictionary.
    """
    design = primer3plus.new()
    design.presets.as_cloning_task()
    if region.direction == -1:
        region = region.flip()
        template = rc(template)

    if lseq and left_overhang:
        raise DasiSequenceDesignException
    if rseq and right_overhang:
        raise DasiSequenceDesignException

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

    if left_overhang is None:
        left_overhang = ""
    if right_overhang is None:
        right_overhang = ""

    design.presets.product_size((len(region), len(region)))
    design.presets.left_overhang(left_overhang)
    design.presets.right_overhang(right_overhang)
    design.PRIMER_PICK_ANYWAY = True
    design.presets.use_overhangs()
    design.presets.long_ok()

    design.logger.set_level("INFO")
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


def no_none_or_nan(*i):
    for _i in i:
        if _i is not None and not np.isnan(_i):
            return _i


def get_primer_extensions(
    graph: nx.DiGraph, n1: AssemblyNode, n2: AssemblyNode, cyclic: bool = True
):
    successors = list(graph.successors(n2))
    if successors:
        sedge = graph[n2][successors[0]]
        r1 = sedge["rprimer_right_ext"]
        r2 = sedge["right_ext"]
        right_ext = no_none_or_nan(r2, r1)
    elif cyclic:
        raise DasiSequenceDesignException
    else:
        right_ext = 0

    predecessors = list(graph.predecessors(n1))
    if predecessors:
        pedge = graph[predecessors[0]][n1]
        l1 = pedge["lprimer_left_ext"]
        l2 = pedge["left_ext"]
        left_ext = no_none_or_nan(l2, l1)
    elif cyclic:
        raise DasiSequenceDesignException
    else:
        left_ext = 0
    return int(left_ext), int(right_ext)


def design_edge(
    assembly: Assembly,
    n1: AssemblyNode,
    n2: AssemblyNode,
    seqdb: Dict[str, SeqRecord],
    query_key: str,
):
    graph = assembly.graph

    edge = n1, n2, graph[n1][n2]

    moltype = edge[2]["type_def"]
    qrecord = seqdb[query_key]
    # contains information about templates and queries

    sequence_result = {}
    if edge[-1]["type_def"].int_or_ext == "external":
        if moltype.use_direct:
            # this is a fragment used directly in an assembly
            sequence = _use_direct(edge, seqdb)
        elif moltype.synthesize:
            # this is either a gene synthesis fragment or already covered by the primers.
            sequence = _design_gap(edge, qrecord)
        else:
            sequence = ""
        sequence_result["sequence"] = sequence
    else:
        pairs, explain = _design_pcr_product_primers(edge, graph, moltype.design, seqdb)
        if pairs:
            sequence_result["primers"] = pairs
        else:
            sequence_result["primers"] = pairs
        sequence_result["primer_explain"] = explain
        if not pairs:
            raise DasiNoPrimerPairsException("No primer pairs were found.")
    return sequence_result


def _use_direct(
    edge: Tuple[AssemblyNode, AssemblyNode, dict], seqdb: Dict[str, SeqRecord]
):
    groups = edge["groups"]
    group = groups[0]
    sk = group.subject_keys[0]
    srecord = seqdb[sk]
    return {
        "SUBJECT_KEY": sk,
        "QUERY_REGION": (group.query_region.a, group.query_region.b),
        "SEQUENCE": str(srecord.seq),
    }


def _skip():
    return {}, {}


def _design_gap(edge: Tuple[AssemblyNode, AssemblyNode, dict], qrecord: SeqRecord):
    n1, _, edge_data = edge
    gene_size = edge_data["gene_size"]
    if not np.isnan(gene_size):
        lshift = edge_data["lshift"]
        assert not np.isnan(lshift)
        a = n1.index + lshift
        b = a + gene_size
        gene_region = edge_data["query_region"].new(a, b)
        gene_seq = gene_region.get_slice(qrecord.seq, as_type=str)
        return {
            "SUBJECT_KEY": None,
            "QUERY_REGION": (gene_region.a, gene_region.b),
            "SEQUENCE": str(gene_seq.seq),
        }
    else:
        return {"SUBJECT_KEY": None, "QUERY_REGION": None, "SEQUENCE": None}


def _design_pcr_product_primers(
    edge: Tuple[AssemblyNode, AssemblyNode, dict],
    graph: nx.DiGraph,
    design: Tuple[bool, bool],
    seqdb: Dict[str, SeqRecord],
):
    if edge[-1]["type_def"].int_or_ext == "external":
        raise Exception()

    # this is a new PCR product
    n1, n2, edata = edge
    lext, rext = get_primer_extensions(graph, n1, n2)
    alignment_groups = edata["groups"]
    group = alignment_groups[0]

    qkey = group.query_key
    qrecord = seqdb[qkey]
    qregion = group.query_region

    # set overhangs
    loverhang = ""
    roverhang = ""
    if design[0]:
        lregion = qregion.new(qregion.a - lext, qregion.a)
        loverhang = lregion.get_slice(qrecord.seq, as_type=str)
    if design[1]:
        rregion = qregion.new(qregion.b, qregion.b + rext)
        roverhang = primer3plus.utils.reverse_complement(
            rregion.get_slice(qrecord.seq, as_type=str)
        )

    # TODO: complex alignment groups have re-adjusted subject regions
    # TODO: ensure fwd and rev primers have same template
    # TODO: the 'fwd' primer of a template in reverse direction is its reverse primer,
    #       so get_template needs to be adjusted.
    # collect template, left primer, and right primer keys

    lkey, rkey = None, None
    if design == (1, 1):
        assert isinstance(group, AlignmentGroup)
        tkey = group.subject_keys[0]
        region = group.alignments[0].subject_region
    else:
        grouping = group.groupings[0]
        template = grouping["template"]
        fwd = grouping["fwd"]
        rev = grouping["rev"]

        tkey = template.subject_key
        if fwd:
            lkey = fwd.subject_key
        if rev:
            rkey = rev.subject_key
        _t = group.get_template(0)
        assert _t.subject_key == tkey
        region = _t.subject_region

    if not design[1]:
        roverhang = ""

    if not design[0]:
        loverhang = ""

    trecord = seqdb[tkey]
    tseq = str(trecord.seq)
    if rkey:
        rrecord = seqdb[rkey]
        rseq = str(rrecord.seq)
    else:
        rseq = None
    if lkey:
        lrecord = seqdb[lkey]
        lseq = str(lrecord.seq)
    else:
        lseq = None

    # design primers
    pairs, explain = design_primers(
        tseq, region, lseq, rseq, left_overhang=loverhang, right_overhang=roverhang
    )
    for pair in pairs:
        pair["LEFT"]["SUBJECT_KEY"] = lkey

        pair["RIGHT"]["SUBJECT_KEY"] = rkey
        pair["PAIR"]["SUBJECT_KEY"] = tkey
    return pairs, explain
