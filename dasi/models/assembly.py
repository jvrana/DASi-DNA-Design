"""Assembly."""
from collections import namedtuple
from copy import deepcopy
from itertools import groupby
from itertools import zip_longest
from typing import Any
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union

import networkx as nx
import numpy as np
import pandas as pd
import primer3plus
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from more_itertools import pairwise
from primer3plus.utils import reverse_complement as rc
from pyblast.utils import is_circular

from .alignment import Alignment
from .alignment import AlignmentGroup
from .alignment import PCRProductAlignmentGroup
from .alignment_container import AlignmentContainer
from .molecule import Molecule
from .molecule import MoleculeType
from .molecule import Reaction
from dasi.config import Config
from dasi.constants import Constants
from dasi.cost import SpanCost
from dasi.exceptions import DasiDesignException
from dasi.exceptions import DasiInvalidMolecularAssembly
from dasi.exceptions import DasiNoPrimerPairsException
from dasi.exceptions import DasiSequenceDesignException
from dasi.log import logger
from dasi.utils import NumpyDataFrame
from dasi.utils import Region
from dasi.utils import sort_cycle


def _design_primers(
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
    design.logger = logger(design)
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
    design.PRIMER_PICK_ANYWAY = False
    design.PRIMER_MIN_ANNEAL_CHECK = Config.PRIMER3_MIN_ANNEAL_CHECK
    design.settings.use_overhangs()
    design.settings.long_ok()

    design.logger.set_level("INFO")

    # TODO: remove debugging code
    try:
        pairs, explain = design.run_and_optimize(
            Config.PRIMER3_N_RUNS, pick_anyway=True
        )
    except Exception as e:
        import json

        print(json.dumps(dict(design.params.items()), indent=2))
        raise e
    if index is not None:
        for pair in pairs.values():
            loc = pair["LEFT"]["location"]
            pair["LEFT"]["location"] = (index[loc[0]], loc[1])

            loc = pair["RIGHT"]["location"]
            pair["RIGHT"]["location"] = (index[loc[0]], loc[1])
    return pairs, explain


def _edata_to_npdf(edata: dict, span_cost: SpanCost) -> NumpyDataFrame:
    return span_cost.cost(edata["span"], edata["type_def"])


def _no_none_or_nan(*i):
    for _i in i:
        if _i is not None and not np.isnan(_i):
            return _i


def _get_primer_extensions(
    graph: nx.DiGraph, n1: "AssemblyNode", n2: "AssemblyNode", cyclic: bool = True
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
    # NOTE: the successors 'left' primer is the fragments 'right' primer extension
    successors = list(graph.successors(n2))
    if successors:
        sedge = graph[n2][successors[0]]
        right_ext = _no_none_or_nan(sedge["lprimer_left_ext"], sedge["left_ext"])
    elif cyclic:
        # TODO: FIX THIS ANNOYING EXCEPTION
        raise DasiSequenceDesignException(
            "Sequence is cyclic but there are no successors for {}".format(n2)
        )
    else:
        right_ext = 0

    # the predecessors 'right' primer is the fragments 'left' primer extension
    predecessors = list(graph.predecessors(n1))
    if predecessors:
        pedge = graph[predecessors[0]][n1]
        left_ext = _no_none_or_nan(pedge["rprimer_right_ext"], pedge["right_ext"])
    elif cyclic:
        # TODO: FIX THIS ANNOYING EXCEPTION
        raise DasiSequenceDesignException(
            "Sequence is cyclic but there are no predecessors for {}".format(n1)
        )
    else:
        left_ext = 0
    return int(left_ext), int(right_ext)


def _use_direct(
    edge: Tuple["AssemblyNode", "AssemblyNode", dict], seqdb: Dict[str, SeqRecord]
) -> Tuple[SeqRecord, AlignmentGroup]:
    group = edge[2]["groups"][0]
    alignment = group.alignments[0]
    sk = alignment.subject_keys
    srecord = seqdb[sk]
    return srecord, alignment


def _design_gap(
    edge: Tuple["AssemblyNode", "AssemblyNode", dict], qrecord: SeqRecord
) -> Union[Tuple[SeqRecord, Region], Tuple[None, None]]:
    n1, _, edge_data = edge
    gene_size = edge_data["gene_size"]
    if not np.isnan(gene_size):
        lshift = edge_data["lshift"]
        assert not np.isnan(lshift)
        a = n1.index + lshift
        b = a + gene_size
        gene_region = edge_data["query_region"].new(a, b)
        gene_seq = gene_region.get_slice(qrecord)
        return gene_seq, gene_region
    else:
        return None, None


# TODO: refactor `_design_pcr_product_primers`, too complex
def _design_pcr_product_primers(
    edge: Tuple["AssemblyNode", "AssemblyNode", dict],
    graph: nx.DiGraph,
    design: Tuple[bool, bool],
    seqdb: Dict[str, SeqRecord],
) -> Tuple[Dict, Dict, AlignmentGroup, Region]:
    if edge[-1]["type_def"].int_or_ext == "external":
        raise Exception()

    # this is a new PCR product
    n1, n2, edata = edge
    lext, rext = _get_primer_extensions(graph, n1, n2)
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
        template = group.get_template()

        # get primer keys
        fwd = group.get_fwd()
        rev = group.get_rev()
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
    try:
        pairs, explain = _design_primers(
            tseq,
            template.subject_region,
            lseq,
            rseq,
            left_overhang=loverhang,
            right_overhang=roverhang,
        )
    except Exception as e:
        raise DasiSequenceDesignException(
            "Could not design primers for {name}[{i}:{j}].\nError: {e}\n{edge}\n{x}".format(
                name=trecord.name,
                i=template.subject_region.a,
                j=template.subject_region.b,
                e=str(e),
                edge=edge,
                x=edge[2]["groups"][0].alignments[0].subject_region,
            )
        ) from e
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


def _design_edge(
    assembly: "Assembly",
    n1: "AssemblyNode",
    n2: "AssemblyNode",
    seqdb: Dict[str, SeqRecord],
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
            return Reaction(
                Reaction.Types.Direct, inputs=[], outputs=[frag_mol], metadata=edge[2]
            )
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
                return Reaction(
                    Reaction.Types.Synthesize,
                    inputs=[],
                    outputs=[synthesis_mol],
                    metadata=edge[2],
                )
            else:
                return None
        else:
            return None
    elif edge[-1]["type_def"].name == Constants.FRAGMENT:
        group = edge[2]["groups"][0]
        alignment = group.alignments[0]
        sk = alignment.subject_key
        frag_seq = seqdb[sk]
        query_region = edge[2]["query_region"]
        frag_mol = Molecule(moltype, alignment, frag_seq, query_region=query_region)
        return Reaction(
            Reaction.Types.Direct, inputs=[], outputs=[frag_mol], metadata=edge[2]
        )
    elif edge[-1]["type_def"].name == Constants.SHARED_SYNTHESIZED_FRAGMENT:
        query_region = edge[2]["query_region"]
        group = edge[2]["groups"][0]
        synthesis_mol = Molecule(
            MoleculeType.types[Constants.SHARED_SYNTHESIZED_FRAGMENT],
            alignment_group=group,
            sequence=query_region.get_slice(seqdb[query_key]),
            query_region=query_region,
        )
        return Reaction(
            Reaction.Types.Synthesize,
            inputs=[],
            outputs=[synthesis_mol],
            metadata=edge[2],
        )
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
            primer_record = SeqRecord(Seq(primer_seq))
            primer = Molecule(
                MoleculeType.types[Constants.PRIMER],
                primer_group,
                primer_record,
                metadata=pair[x],
            )
            primers.append(primer)
        template = Molecule(
            MoleculeType.types[Constants.TEMPLATE],
            pair["PAIR"]["GROUP"],
            seqdb[pair["PAIR"]["SUBJECT_KEY"]],
        )

        product = Molecule(
            moltype,
            group,
            query_region.get_slice(seqdb[query_key]),
            query_region=query_region,
        )

        return Reaction(
            Reaction.Types.PCR,
            inputs=primers + [template],
            outputs=[product],
            metadata=edge[2],
        )


AssemblyNode = namedtuple(
    "AssemblyNode", "index expandable type overhang"
)  #: tuple representing a location on a goal sequence


class Assembly(Iterable):
    """Should take in a path, graph, container, seqdb to produce relevant
    information."""

    REACTION_DF_COLS = (
        "DESIGN_ID",
        "DESIGN_KEY",
        "ASSEMBLY_ID",
        "REACTION_ID",
        "REACTION_NAME",
        "NAME",
        "TYPE",
        "KEY",
        "ROLE",
        "REGION",
        "SEQUENCE",
        "LENGTH",
        "META",
    )
    REACTION_DF_SORT_BY = (
        "TYPE",
        "DESIGN_ID",
        "ASSEMBLY_ID",
        "REACTION_ID",
        "REACTION_NAME",
        "NAME",
        "ROLE",
    )
    SUMMARY_DF_COLS = (
        "query_start",
        "query_end",
        "subject_names",
        "subject_keys",
        "subject_starts",
        "subject_ends",
        "cost",
        "material",
        "span",
        "type",
        "internal_or_external",
        "efficiency",
        "complexity",
        "notes",
    )

    def __init__(
        self,
        nodes: List[AssemblyNode],
        container: AlignmentContainer,
        full_assembly_graph: nx.DiGraph,
        query_key: str,
        query: SeqRecord,
        seqdb,
        do_raise: bool = True,
    ):
        self.logger = logger(self)
        self._nodes = tuple(nodes)
        self._reactions = tuple()
        self.validate_input_nodes()
        self.seqdb = seqdb
        self.query_key = query_key
        self.query = query
        self._full_graph = full_assembly_graph

        self.container = container
        self.groups = container.groups()
        if len(self.groups) == 0:
            raise DasiDesignException("No groups were found in container.")
        self.graph = self._subgraph(self._full_graph, nodes, do_raise=do_raise)
        nx.freeze(self.graph)

        if do_raise:
            self.post_validate(do_raise)

    @property
    def reactions(self):

        top_level_reaction = Reaction(
            Reaction.Types.Assembly,
            inputs=[],
            outputs=[],
            metadata={"query_key": self.query_key},
        )
        top_level_reaction.outputs = [
            Molecule(
                molecule_type=MoleculeType.types[Constants.PLASMID],
                sequence=deepcopy(self.query),
                alignment_group=None,
            )
        ]

        if not self._reactions:
            reactions = [top_level_reaction]
            for n1, n2, edata in self.edges():
                reaction = _design_edge(self, n1, n2, seqdb=self.seqdb)
                if reaction:
                    for _mol in reaction.outputs:
                        top_level_reaction.inputs.append(_mol)
                    reactions.append(reaction)
            self._reactions = tuple(reactions)
        return self._reactions

    @property
    def assembly_reactions(self):
        return [r for r in self.reactions if r.name == Reaction.Types.Assembly]

    @property
    def nonassembly_reactions(self):
        return [r for r in self.reactions if r.name != Reaction.Types.Assembly]

    def post_validate(self):
        total_span = 0
        for n1, n2, edata in self.edges():
            if n1.type == n2.type:
                raise DasiInvalidMolecularAssembly("Invalid assembly graph")
            total_span += edata["span"]
        if not total_span == len(self.query):
            raise DasiInvalidMolecularAssembly(
                "Assembly length '{}' is different from expected"
                " length '{}'".format(total_span, len(self.query))
            )

    def _head(self):
        """Get the 'first' 'A' node."""
        x = sorted(list(self.graph.nodes), key=lambda n: (n.type == "B", n.index))
        return x[0]

    def _sorted_edges(self):
        head = self._head()
        edges = list(nx.bfs_edges(self.graph, head))
        edges += [
            (t[1], t[0])
            for t in nx.bfs_edges(self.graph, head, reverse=True, depth_limit=1)
        ]

        return edges

    def validate_input_nodes(self):
        # rule 1 A -> B -> A -> B -> ...
        types = [n.type for n in self._nodes]
        groups = [list(g) for _, g in groupby(types)]
        if len(types) != len(groups):
            raise ValueError("There invalid edges input nodes")

    @staticmethod
    def _missing_edata():
        return {
            "cost": np.inf,
            "weight": np.inf,
            "material": np.inf,
            "efficiency": 0.0,
            "type_def": MoleculeType.types[Constants.MISSING],
            "span": np.inf,
        }

    def _subgraph(
        self, graph: nx.DiGraph, nodes: List[AssemblyNode], do_raise: bool = True
    ):
        def _resolve(node: "AssemblyNode", qregion) -> Tuple[AssemblyNode, dict]:
            new_node = AssemblyNode(qregion.t(node.index), *list(node)[1:])
            return new_node

        subgraph = nx.OrderedDiGraph()
        nodes = [AssemblyNode(*n) for n in nodes]
        example_query_region = self.container.alignments[0].query_region

        resolved_nodes = [_resolve(node, example_query_region) for node in nodes]
        if self.cyclic:
            resolved_nodes = sort_cycle(
                resolved_nodes, key=lambda n: (n.type, n.index, n)
            )
        subgraph.add_nodes_from(resolved_nodes)

        pair_iter = list(pairwise(nodes))

        if self.cyclic:
            pair_iter.append((nodes[-1], nodes[0]))

        visited = set()
        for n1, n2 in pair_iter:
            edata = graph.get_edge_data(n1, n2)
            if edata is None:
                edata = self._missing_edata()
            else:
                assert edata["type_def"].int_or_ext

            query_region = self.container.alignments[0].query_region.new(
                n1.index, n2.index
            )
            groups = self.container.find_groups_by_pos(
                query_region.a,
                query_region.b,
                group_type=edata["type_def"].name,
                groups=self.groups,
            )
            if edata["type_def"].int_or_ext == "internal" and not groups:
                raise DasiDesignException(
                    "Missing groups for edge between {} and {}".format(n1, n2)
                )

            edata["groups"] = groups
            edata["query_region"] = query_region

            rn1 = _resolve(n1, query_region)
            rn2 = _resolve(n2, query_region)

            if rn1 in visited:
                break
            visited.add(rn1)

            # TODO: add this check
            if do_raise:
                if rn1 in subgraph:
                    raise DasiInvalidMolecularAssembly(
                        "Node already exists in subgraph"
                    )
                if rn2 in subgraph:
                    raise DasiInvalidMolecularAssembly(
                        "Node already exists in subgraph"
                    )

            subgraph.add_edge(rn1, rn2, **edata)
        return subgraph

    @property
    def cyclic(self):
        return is_circular(self.query)

    # TODO: consolidate this with shortest path utils in networkx
    def compute_cost(self):
        c = self._compute_cost()
        if np.isinf(c["material"]):
            return np.inf
        return c["material"] / c["efficiency"]

    def _compute_cost(self):
        material = 0
        efficiency = 1.0
        for _, _, edata in self.edges():
            material += edata["material"]
            efficiency *= edata["efficiency"]
        if efficiency == 0:
            return {"material": np.inf, "efficiency": 0.0}
        return {"material": material, "efficiency": efficiency}

    @property
    def cost(self):
        return self.compute_cost()

    @property
    def material_cost(self):
        return self._compute_cost()["material"]

    @property
    def efficiency(self):
        return self._compute_cost()["efficiency"]

    def edges(self, data=True) -> Iterable[Tuple[AssemblyNode, AssemblyNode, Dict]]:
        for n1, n2 in self._sorted_edges():
            if data:
                edata = self.graph[n1][n2]
                yield n1, n2, edata
            else:
                yield n1, n2

    def nodes(self, data=True) -> Iterable[Tuple[AssemblyNode, Dict]]:
        return self.graph.nodes(data=data)

    def edit_distance(
        self, other: "Assembly", explain=False
    ) -> Union[int, Tuple[int, List[Tuple[int, str]]]]:
        differences = []
        for i, (n1, n2) in enumerate(
            zip_longest(self.nodes(data=False), other.nodes(data=False))
        ):
            if n1 is None or n2 is None:
                differences.append((i, "{} != {}".format(n1, n2)))
                continue
            if n1.index != n2.index:
                differences.append((i, "Index: {} != {}".format(n1.index, n2.index)))
            if n1.expandable != n2.expandable:
                differences.append(
                    (i, "Expandable: {} != {}".format(n1.expandable, n2.expandable))
                )
            if n1.type != n2.type:
                differences.append((i, "Type: {} != {}".format(n1.type, n2.type)))
            if n1.overhang != n2.overhang:
                differences.append(
                    (i, "Overhang: {} != {}".format(n1.overhang, n2.overhang))
                )
        dist = len(differences)
        if explain:
            return dist, differences
        return dist

    def print(self):
        print("query_name: {}".format(self.query.name))
        print("query_key: {}".format(self.query_key))
        print("Cost: {}".format(self.compute_cost()))
        df = self.to_df()
        print(df)

    def print_diff(self, other: "Assembly"):
        for i, (n1, n2) in enumerate(
            zip_longest(self.nodes(data=False), other.nodes(data=False))
        ):
            if n1 != n2:
                desc = False
            else:
                desc = True
            print("{} {} {}".format(desc, n1, n2))

    def to_df(self) -> pd.DataFrame:
        rows = []
        for n1, n2, edata in self.edges():
            groups = edata["groups"]

            if groups:
                group = groups[0]
                if isinstance(group, PCRProductAlignmentGroup):
                    alignments = group.alignments
                else:
                    alignments = group.alignments[:1]
            else:
                alignments = []
            subject_keys = [a.subject_key for a in alignments]
            subject_names = [self.container.seqdb[key].name for key in subject_keys]
            subject_starts = [a.subject_region.a for a in alignments]
            subject_ends = [a.subject_region.b for a in alignments]

            rows.append(
                {
                    "query_start": edata["query_region"].a,
                    "query_end": edata["query_region"].b,
                    "subject_names": subject_names,
                    "subject_keys": subject_keys,
                    "subject_starts": subject_starts,
                    "subject_ends": subject_ends,
                    "cost": edata["cost"],
                    "material": edata["material"],
                    "span": edata["span"],
                    "type": edata["type_def"].name,
                    "internal_or_external": edata["type_def"].int_or_ext,
                    "efficiency": edata.get("efficiency", np.nan),
                    "complexity": edata.get("complexity", np.nan),
                    "notes": edata.get("notes", ""),
                }
            )

        df = pd.DataFrame(rows)
        assert set(df.columns) == set(self.SUMMARY_DF_COLS)
        return df

    def _csv_row(
        self,
        m: Molecule,
        role: str,
        reaction_id: Union[str, int],
        reaction_name: str,
        meta: dict = None,
    ) -> dict:
        mtype = m.type.name
        group = m.alignment_group

        key = None
        name = None
        if group:
            if hasattr(group, "subject_key"):
                if group.subject_key:
                    key = group.subject_key
                    name = self.seqdb[key].name
            else:
                # TODO: select best subject here
                key = group.subject_keys[0]
                name = self.seqdb[key].name

        if m.query_region:
            q = (m.query_region.a, m.query_region.b, m.query_region.context_length)
        else:
            q = None
        length = len(m.sequence)
        data = {
            "NAME": name,
            "LENGTH": length,
            "SEQUENCE": str(m.sequence.seq),
            "TYPE": mtype,
            "KEY": key,
            "REGION": q,
            "ROLE": role,
            "REACTION_ID": reaction_id,
            "REACTION_NAME": reaction_name,
        }
        if meta:
            data.update({"META": deepcopy(meta)})
        return data

    @property
    def molecules(self) -> Generator[Tuple[int, str, Molecule], None, None]:
        for i, r in enumerate(self.reactions):
            for m in r.inputs:
                yield (i, "input", m)
            for m in r.outputs:
                yield (i, "output", m)

    def to_reaction_df(self) -> pd.DataFrame:
        rows = []
        reactions = self.reactions
        for reaction_id, role, m in self.molecules:
            reaction = reactions[reaction_id]
            if m.type.name == "PRIMER":
                meta = deepcopy(m.metadata)
                meta["ANNEAL"] = meta["SEQUENCE"]
                del meta["SEQUENCE"]
                meta = {"PRIMER_{}".format(k): v for k, v in meta.items()}
            else:
                meta = None
            rows.append(self._csv_row(m, role, reaction_id, reaction.name, meta))
        colnames = self.REACTION_DF_COLS
        df = pd.DataFrame(rows, columns=colnames)
        df.sort_values(by=list(self.REACTION_DF_SORT_BY), inplace=True)
        return df

    def __eq__(self, other: "Assembly") -> bool:
        return self.edit_distance(other) == 0

    def __iter__(self) -> Generator[AssemblyNode, None, None]:
        yield from self.nodes(data=False)
