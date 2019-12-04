from typing import Any
from typing import Dict
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
import primer3plus
from Bio.SeqRecord import SeqRecord
from primer3plus.utils import reverse_complement as rc

from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.exceptions import DasiNoPrimerPairsException
from dasi.exceptions import DasiSequenceDesignException
from dasi.models import Alignment
from dasi.models import AlignmentGroup
from dasi.models import Assembly
from dasi.models import AssemblyNode
from dasi.models import Molecule
from dasi.models import MoleculeType
from dasi.models import Reaction
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
    :param left_overhang: optionally provided left overhang sequence of the primer
    :param right_overhang: optionally provided right overhang sequence of the primer
    :return: tuple of pairs and the 'explain' dictionary.
    """
    design = primer3plus.new()
    design.settings.as_cloning_task()
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
        design.settings.template(adjusted_template)
        design.settings.included((0, len(region)))
        index = list(region) + list(region.invert()[0])
    else:
        design.settings.template(template)
        design.settings.included((region.a, len(region)))
        index = None
    if lseq:
        design.settings.left_sequence(lseq)
    if rseq:
        design.settings.right_sequence(rseq)

    if left_overhang is None:
        left_overhang = ""
    if right_overhang is None:
        right_overhang = ""

    design.settings.product_size((len(region), len(region)))
    design.settings.left_overhang(left_overhang)
    design.settings.right_overhang(right_overhang)
    design.PRIMER_PICK_ANYWAY = True
    design.settings.use_overhangs()
    design.settings.long_ok()

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
) -> Tuple[int, int]:
    """Return the left and right primer extensions for the given *internal*
    fragment. To get the extensions, we look for the left predecessor edge and
    get its `right_ext` or `rprimer_right_ext` and on the other side the right
    successor edge and its `left_ext` or `lprimer_left_ext`.

    :param graph: assembly graph
    :param n1: source node on the graph
    :param n2: end node on the graph
    :param cyclic: whether the
    :return: tuple of left and right extensions
    """
    # the successors 'left' primer is the fragments 'right' primer extension
    successors = list(graph.successors(n2))
    if successors:
        sedge = graph[n2][successors[0]]
        right_ext = no_none_or_nan(sedge["lprimer_left_ext"], sedge["left_ext"])
    elif cyclic:
        raise DasiSequenceDesignException
    else:
        right_ext = 0

    # the predecessors 'right' primer is the fragments 'left' primer extension
    predecessors = list(graph.predecessors(n1))
    if predecessors:
        pedge = graph[predecessors[0]][n1]
        left_ext = no_none_or_nan(pedge["rprimer_right_ext"], pedge["right_ext"])
    elif cyclic:
        raise DasiSequenceDesignException
    else:
        left_ext = 0
    return int(left_ext), int(right_ext)


def _use_direct(
    edge: Tuple[AssemblyNode, AssemblyNode, dict], seqdb: Dict[str, SeqRecord]
) -> Tuple[str, AlignmentGroup]:
    group = edge[2]["groups"][0]
    alignment = group.alignments[0]
    sk = alignment.subject_keys
    srecord = seqdb[sk]
    return str(srecord.seq), alignment


def _design_gap(
    edge: Tuple[AssemblyNode, AssemblyNode, dict], qrecord: SeqRecord
) -> Union[Tuple[str, Region], Tuple[None, None]]:
    n1, _, edge_data = edge
    gene_size = edge_data["gene_size"]
    if not np.isnan(gene_size):
        lshift = edge_data["lshift"]
        assert not np.isnan(lshift)
        a = n1.index + lshift
        b = a + gene_size
        gene_region = edge_data["query_region"].new(a, b)
        gene_seq = gene_region.get_slice(qrecord.seq)
        return gene_seq, gene_region
    else:
        return None, None


def _design_pcr_product_primers(
    edge: Tuple[AssemblyNode, AssemblyNode, dict],
    graph: nx.DiGraph,
    design: Tuple[bool, bool],
    seqdb: Dict[str, SeqRecord],
) -> Tuple[Dict, Dict, AlignmentGroup, Region]:
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
    fwd, template, rev = None, None, None
    if design == (1, 1):
        assert isinstance(group, AlignmentGroup)
        tkey = group.subject_keys[0]
        template = group.alignments[0]
    else:
        grouping = group.groupings[0]
        template = group.get_template(0)
        assert grouping["template"].subject_key == template.subject_key

        # get primer keys
        fwd = grouping["fwd"]
        rev = grouping["rev"]
        tkey = template.subject_key
        if fwd:
            lkey = fwd.subject_key
        if rev:
            rkey = rev.subject_key

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
        tseq,
        template.subject_region,
        lseq,
        rseq,
        left_overhang=loverhang,
        right_overhang=roverhang,
    )
    for pair in pairs.values():
        pair["LEFT"]["SUBJECT_KEY"] = lkey
        pair["LEFT"]["GROUP"] = fwd
        pair["RIGHT"]["SUBJECT_KEY"] = rkey
        pair["RIGHT"]["GROUP"] = rev
        pair["PAIR"]["SUBJECT_KEY"] = tkey
        pair["PAIR"]["GROUP"] = template

    query_region = group.query_region.new(
        group.query_region.a - len(loverhang) + group.query_region.context_length,
        group.query_region.c + len(roverhang) + group.query_region.context_length,
    )
    #     print(group.query_region)
    #     print(template.query_region)
    #     print(group)
    #     print(len(query_region), len(group.query_region) + lext + rext, lext, rext)
    #     assert len(query_region) == len(group.query_region) + lext + rext
    return pairs, explain, template, query_region


def design_edge(
    assembly: Assembly, n1: AssemblyNode, n2: AssemblyNode, seqdb: Dict[str, SeqRecord]
) -> Union[Reaction, None]:
    query_key = assembly.query_key
    graph = assembly.graph

    edge = n1, n2, graph[n1][n2]

    moltype = edge[2]["type_def"]
    qrecord = seqdb[query_key]
    # contains information about templates and queries

    if edge[-1]["type_def"].int_or_ext == "external":
        if moltype.use_direct:
            # this is a fragment used directly in an assembly
            frag_seq, frag_alignment = _use_direct(edge, seqdb)
            frag_mol = Molecule(moltype, frag_alignment, frag_seq)
            return Reaction("Use Direct", inputs=[], outputs=[frag_mol])
        elif moltype.synthesize:
            # this is either a gene synthesis fragment or already covered by the primers.
            synthesis_seq, synthesis_region = _design_gap(edge, qrecord)
            if synthesis_seq:
                subject_region = Region(
                    0, len(synthesis_region), len(synthesis_region), direction=1
                )
                synthesis_alignment = Alignment(
                    synthesis_region, subject_region, moltype.name, query_key, ""
                )
                synthesis_mol = Molecule(
                    moltype,
                    synthesis_alignment,
                    synthesis_seq,
                    query_region=synthesis_region,
                )
                return Reaction("Synthesize", inputs=[], outputs=[synthesis_mol])
            else:
                return None
        else:
            return None
    else:
        pairs, explain, group, query_region = _design_pcr_product_primers(
            edge, graph, moltype.design, seqdb
        )
        if not pairs:
            raise DasiNoPrimerPairsException("No primer pairs were found.")
        pair = pairs[0]
        primers = []
        for x in ["LEFT", "RIGHT"]:
            primer_seq = pair[x]["OVERHANG"] + pair[x]["SEQUENCE"]
            primer_group = pair[x]["GROUP"]
            primer = Molecule(
                MoleculeType.types[Constants.PRIMER],
                primer_group,
                primer_seq,
                metadata=pair[x],
            )
            primers.append(primer)
        template = Molecule(
            MoleculeType.types[Constants.TEMPLATE],
            pair["PAIR"]["GROUP"],
            str(seqdb[pair["PAIR"]["SUBJECT_KEY"]].seq),
        )

        product = Molecule(
            moltype,
            group,
            query_region.get_slice(seqdb[query_key]),
            query_region=query_region,
        )

        return Reaction("PCR", inputs=primers + [template], outputs=[product])
