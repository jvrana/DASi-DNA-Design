import hashlib
from typing import Dict
from typing import Tuple
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


def dasi_design_to_dag(design: Union["Design", "LibraryDesign"]) -> nx.DiGraph:
    # TODO: standard way to display SeqRecord
    graph = nx.DiGraph()
    for q_i, (qk, result) in enumerate(design.results.items()):
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
    return graph


def _validate_dag_graph(graph: nx.DiGraph):
    """Validates the DASi DAG graph."""
    for n2, ndata2 in graph.nodes(data=True):
        ########################
        # validate each non-plasmid is expected to
        # have at most one reaction
        ########################
        molecule = ndata2.get("molecule", None)
        if molecule and molecule.type is not MoleculeType.types[C.PLASMID]:
            predecessors = list(graph.predecessors(n2))
            if len(predecessors) > 1:
                msg = str(n2)
                for n1 in graph.predecessors(n2):
                    msg += "\n\t" + str(n1)
                raise ValueError(
                    "Molecule {} has more than one reaction.\n{}".format(n2, msg)
                )

        ########################
        # validate each reaction
        # has expected number of
        # predecessors and successors
        ########################
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


def _used_in_designs(ndata: Dict) -> Dict:
    x = sorted(list(ndata["design_keys"]))
    return [{"design_key": _x[0], "assembly": _x[1]} for _x in x]


def _reactions_property(
    graph: nx.DiGraph,
    reaction_node_dict: Dict[str, int],
    molecule_node_dict: Dict[str, int],
) -> Dict:
    property_reaction = []
    for n, index in reaction_node_dict.items():
        ndata = graph.nodes[n]
        output_keys = sorted([mhash(m) for m in ndata["reaction"].outputs])
        input_keys = sorted([mhash(m) for m in ndata["reaction"].inputs])
        property_reaction.append(
            {
                "name": ndata["reaction"].name,
                "index": index,
                "used_in_assemblies": _used_in_designs(ndata),
                "inputs": [molecule_node_dict[n] for n in input_keys],
                "outputs": [molecule_node_dict[n] for n in output_keys],
            }
        )
    return property_reaction


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
                "name": mol.type.name,
                "index": index,
                "sequence": seqrecord_to_json(mol.sequence),
                "used_in_assemblies": _used_in_designs(ndata),
                "used_as_input_to_reactions": used_as_inputs,
                "used_as_output_to_reactions": used_as_outputs,
            }
        )
    return property_molecule


def _design_property(design, reaction_node_dict):
    status = design.status

    for qk, qstatus in status.items():
        qstatus["sequence"] = seqrecord_to_json(design.seqdb[qk])
        for a_i, adata in enumerate(qstatus["assemblies"]):
            assembly = design.results[qk].assemblies[a_i]
            adata["nonassembly_reactions"] = [
                reaction_node_dict[rhash(r, a_i)]
                for r in assembly.nonassembly_reactions
            ]
            adata["assembly_reactions"] = [
                reaction_node_dict[rhash(r, a_i)] for r in assembly.assembly_reactions
            ]

            # TODO: run start
    return status


def dasi_design_to_output_json(design: Union["Design", "LibraryDesign"]):
    graph = dasi_design_to_dag(design)
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

    return {
        "metadata": design.metadata,
        "designs": _design_property(design, reaction_node_dict),
        "molecules": _molecules_property(graph, reaction_node_dict, molecule_node_dict),
        "reactions": _reactions_property(graph, reaction_node_dict, molecule_node_dict),
    }
