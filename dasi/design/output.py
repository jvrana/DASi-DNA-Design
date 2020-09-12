import functools
import hashlib
import json
import operator
from copy import deepcopy
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import TypeVar
from typing import Union

import networkx as nx

from dasi.constants import Constants as C
from dasi.exceptions import DasiOutputValidationError
from dasi.models import Molecule
from dasi.models import MoleculeType
from dasi.models import Reaction
from dasi.schemas import Schemas
from dasi.schemas import validate_with_schema
from dasi.utils.biopython import seqrecord_to_json

DesignType = TypeVar("Design")


# TODO: Add "Order Primer" to list of reactions
def validate_output(data, do_raise: bool = True):
    return validate_with_schema(
        data,
        Schemas.output_schema,
        do_raise=do_raise,
        reraise_as=DasiOutputValidationError,
    )


def shash(seq: str) -> str:
    """Convert sequence string into a hash."""
    return hashlib.sha1(seq.strip().upper().encode()).hexdigest()


def mhash(mol1: Molecule) -> Tuple:
    """Convert molecule to hash."""
    return ("molecule", mol1.type.name, shash(str(mol1.sequence.seq)))


def rhash(rxn: Reaction, assembly_id: int) -> Tuple:
    """Convert reaction to hash."""
    outputs = sorted([mhash(m) for m in rxn.outputs])
    inputs = sorted([mhash(m) for m in rxn.inputs])
    if rxn.name == "Assembly":
        x = ["reaction", rxn.name, assembly_id]
    else:
        x = ["reaction", rxn.name]
    return tuple(x + [tuple(inputs)] + [tuple(outputs)])


def dasi_design_to_dag(
    design: DesignType,
    validate: bool = True,
    elim_extra_reactions: bool = False,
    query_keys: Optional[List[str]] = None,
) -> nx.DiGraph:
    # TODO: standard way to display SeqRecord
    graph = nx.DiGraph()
    for q_i, (qk, result) in enumerate(design.results.items()):
        if query_keys and qk not in query_keys:
            continue
        for a_i, a in enumerate(result.assemblies):
            dk = (qk, a_i)
            for r in a.reactions:
                r1 = rhash(r, a_i)
                graph.add_node(r1, reaction=r, assembly_index=a_i, assembly=a)
                graph.nodes[r1].setdefault("design_keys", set()).add(dk)
                for in_mol in r.inputs:
                    n1 = mhash(in_mol)
                    graph.add_node(n1, molecule=in_mol)
                    graph.nodes[n1].setdefault("design_keys", set()).add(dk)
                    graph.add_edge(n1, r1)
                for out_mol in r.outputs:
                    n2 = mhash(out_mol)
                    graph.add_node(n2, molecule=out_mol)
                    graph.nodes[n2].setdefault("design_keys", set()).add(dk)
                    graph.add_edge(r1, n2)
    if validate:
        _validate_dag_graph(graph, elim_extra_reactions=elim_extra_reactions)
    return graph


def _validate_num_reactions(graph: nx.DiGraph, elim_extra_reactions: bool = False):
    """Validate that each DNA fragment for each plasmid has exactly one
    reaction. Optionally, remove these reactions from the final DAG.

    :param graph:
    :param elim_extra_reactions:
    :return:
    """
    nodes_to_remove = []
    for n2, ndata2 in graph.nodes(data=True):
        ########################
        # validate each non-plasmid is expected to
        # have at most one reaction
        ########################
        molecule = ndata2.get("molecule", None)

        if molecule and molecule.type is not MoleculeType.types[C.PLASMID]:
            predecessors = list(graph.predecessors(n2))
            if len(predecessors) > 1:
                if elim_extra_reactions:
                    for n1 in list(graph.predecessors(n2))[1:]:
                        nodes_to_remove.append(n1)
                else:
                    msg = str(n2)
                    for n1 in graph.predecessors(n2):
                        msg += "\n\t" + str(n1)
                    raise ValueError(
                        "Molecule {} has more than one reaction.\n{}".format(n2, msg)
                    )

    if nodes_to_remove:
        for n1 in nodes_to_remove:
            graph.remove_node(n1)


def _validate_reactions(graph):
    """Validate each reaction has expected number of predecessors and
    successors."""
    for n2, ndata2 in graph.nodes(data=True):
        reaction = ndata2.get("reaction", None)
        if reaction:
            predecessors = list(graph.predecessors(n2))
            if not len(predecessors) == len(reaction.inputs):
                raise ValueError(
                    "The number of predecessors ({}) does not"
                    " match the number of inputs for the reaction"
                    " ({})".format(len(predecessors), len(reaction.outputs))
                )

            successors = list(graph.successors(n2))
            if not len(successors) == len(reaction.outputs):
                raise ValueError(
                    "The number of successors ({}) does not"
                    " match the number of outputs for the reaction"
                    " ({})".format(len(successors), len(reaction.outputs))
                )


def _validate_dag_graph(graph: nx.DiGraph, elim_extra_reactions: bool = False):
    """Validates the DASi DAG graph."""
    _validate_num_reactions(graph, elim_extra_reactions)
    _validate_reactions(graph)


def _used_in_designs(ndata: Dict) -> Dict:
    x = sorted(list(ndata["design_keys"]))
    return [{"design_key": _x[0], "assembly": _x[1]} for _x in x]


def _reaction_metadata(reaction):
    return {
        "material_cost": reaction.metadata["material"],
        "efficiency": reaction.metadata.get("efficiency", 0.0),
        "complexity": reaction.metadata.get("complexity", None),
    }


def _assembly_metadata(assembly):
    rmetas = [_reaction_metadata(r) for r in assembly.nonassembly_reactions]
    material = sum([m["material_cost"] for m in rmetas])
    efficiency = functools.reduce(operator.mul, [m["efficiency"] for m in rmetas])
    complexities = [m["complexity"] for m in rmetas if m["complexity"] is not None]

    if complexities:
        max_complexity = max(complexities)
    else:
        max_complexity = 0

    return {
        "material_cost": material,
        "assembly_efficiency": efficiency,
        "max_complexity": max_complexity,
    }


def _reactions_property(
    graph: nx.DiGraph,
    reaction_node_dict: Dict[str, int],
    molecule_node_dict: Dict[str, int],
) -> Dict:
    property_reaction = []
    for n, index in reaction_node_dict.items():
        ndata = graph.nodes[n]
        reaction = ndata["reaction"]
        output_keys = sorted([mhash(m) for m in reaction.outputs])
        input_keys = sorted([mhash(m) for m in reaction.inputs])
        if reaction.name != Reaction.Types.Assembly:
            metadata = _reaction_metadata(reaction)
        else:
            assembly = ndata["assembly"]
            metadata = _assembly_metadata(assembly)

        property_reaction.append(
            {
                "__name__": ndata["reaction"].name,
                "__index__": index,
                "__type__": "reaction",
                "used_in_assemblies": _used_in_designs(ndata),
                "inputs": [molecule_node_dict[n] for n in input_keys],
                "outputs": [molecule_node_dict[n] for n in output_keys],
                "metadata": metadata,
            }
        )
    return property_reaction


def _clean_metadata(metadata):
    cleaned = {}
    for k, v in metadata.items():
        if isinstance(v, dict):
            cleaned[k] = _clean_metadata(v)
        else:
            try:
                json.dumps(v)
                cleaned[k] = v
            except Exception:
                pass
    return cleaned


def _molecules_property(
    graph: nx.DiGraph,
    reaction_node_dict: Dict[str, int],
    molecule_node_dict: Dict[str, int],
) -> Dict:
    property_molecule = []
    for n, index in molecule_node_dict.items():
        ndata = graph.nodes[n]
        mol = ndata["molecule"]
        used_as_inputs = [reaction_node_dict[n] for n in graph.successors(n)]
        used_as_outputs = [reaction_node_dict[n] for n in graph.predecessors(n)]
        property_molecule.append(
            {
                "__name__": mol.type.name,
                "__index__": index,
                "__type__": "molecule",
                "__meta__": deepcopy(_clean_metadata(mol.metadata)),
                "sequence": seqrecord_to_json(mol.sequence),
                "used_in_assemblies": _used_in_designs(ndata),
                "used_as_input_to_reactions": used_as_inputs,
                "used_as_output_to_reactions": used_as_outputs,
            }
        )
    return property_molecule


def _design_property(design, reaction_node_dict, graph):
    status = design.status

    def _reaction_summ(r, a_i):
        return {
            "reaction_index": reaction_node_dict.get(rhash(r, a_i)),
            "metadata": {
                "cost": r.metadata["cost"],
                "materials": r.metadata["cost"],
                "efficiency": r.metadata.get("efficiency", 0.0),
                "complexity": r.metadata.get("complexity", None),
            },
            "outputs": [
                {
                    "start": m.query_region.a,
                    "end": m.query_region.b,
                    "type": m.type.name,
                }
                for m in r.outputs
            ],
        }

    for qk, qstatus in status.items():
        qstatus["sequence"] = seqrecord_to_json(design.seqdb[qk])
        for a_i, adata in enumerate(qstatus["assemblies"]):
            assembly = design.results[qk].assemblies[a_i]

            nonassembly_reactions = [
                r for r in assembly.nonassembly_reactions if rhash(r, a_i) in graph
            ]
            assembly_reactions = [
                r for r in assembly.assembly_reactions if rhash(r, a_i) in graph
            ]

            adata["summary"] = [_reaction_summ(r, a_i) for r in nonassembly_reactions]
            adata["final_assembly_reaction"] = [
                reaction_node_dict[rhash(r, a_i)] for r in assembly_reactions
            ]

            # TODO: run start
    return status


def dasi_design_to_output_json(
    design: DesignType,
    elim_extra_reactions: bool = False,
    query_keys: Optional[List[str]] = None,
):
    """Convert a DASi Design instance into an output JSON."""
    graph = dasi_design_to_dag(
        design, elim_extra_reactions=elim_extra_reactions, query_keys=query_keys
    )
    reaction_node_dict = {}
    molecule_node_dict = {}
    sorted_nodes = list(nx.topological_sort(graph))[::-1]
    for n in sorted_nodes:
        ndata = graph.nodes[n]
        r = ndata.get("reaction", None)
        m = ndata.get("molecule", None)
        if r:
            ndata["index"] = len(reaction_node_dict)
            assert n not in reaction_node_dict
            reaction_node_dict[n] = ndata["index"]
        elif m:
            ndata["index"] = len(molecule_node_dict)
            assert n not in molecule_node_dict
            molecule_node_dict[n] = ndata["index"]

            if m.type.name == C.PLASMID:
                pass
        else:
            raise ValueError

    output = {
        "metadata": design.metadata,
        "designs": _design_property(design, reaction_node_dict, graph),
        "molecules": _molecules_property(graph, reaction_node_dict, molecule_node_dict),
        "reactions": _reactions_property(graph, reaction_node_dict, molecule_node_dict),
    }

    OutputValidator.validate_design_assemblies(output)
    return output


def _filter_dict(data, keys):
    return {k: deepcopy(v) for k, v in data.items() if k in keys}


def filter_output_by_design_keys(
    output, keys, only_assemblies: Optional[List[str]] = None
):
    new_out = {}

    new_out["metadata"] = deepcopy(output["metadata"])
    new_out["designs"] = _filter_dict(output["designs"], keys)

    rids = []
    for design in new_out["designs"].values():
        for assembly in design["assemblies"]:
            for reaction in assembly["summary"]:
                rids.append(reaction["reaction_index"])
            rids += assembly["final_assembly_reaction"]
    rids = sorted(list(set(rids)))

    rdict = {r["__index__"]: r for r in output["reactions"]}
    mdict = {r["__index__"]: r for r in output["molecules"]}

    reactions = [rdict[r] for r in rids]

    mids = []
    for r in reactions:
        mids += r["inputs"]
        mids += r["outputs"]
    mids = sorted(list(set(mids)))

    molecules = [mdict[m] for m in mids]

    new_out["reactions"] = reactions
    new_out["molecules"] = molecules
    validate_output(new_out)
    return new_out


class OutputValidator:
    @classmethod
    def validate_design_assemblies(cls, out):
        for d in out["designs"].values():
            for a in d["assemblies"]:
                for rid in a["final_assembly_reaction"]:
                    reaction = out["reactions"][rid]
                    assert reaction["__name__"] == "Assembly"

    # TODO: In Silico from output json, compare to original sequences
