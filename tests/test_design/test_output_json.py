import hashlib
import json
from itertools import product

import networkx as nx

from dasi import Design
from dasi.constants import Constants as C
from dasi.models import MoleculeType
from dasi.models import Reaction


def test_output():
    design = Design.fake(
        n_designs=3, n_linear_seqs=500, n_primers_from_templates=50, shared_length=0
    )
    design.run()

    def seq_sha1(seq: str) -> str:
        """Convert sequence string into a hash."""
        return hashlib.sha1(seq.strip().upper().encode()).hexdigest()

    def mol_key(mol1):
        return ("molecule", mol1.type.name, seq_sha1(str(mol1.sequence.seq)))

    def rxn_key(rxn, assembly_id):
        outputs = sorted([mol_key(m) for m in rxn.outputs])
        inputs = sorted([mol_key(m) for m in rxn.inputs])
        if rxn.name == "Assembly":
            x = ["reaction", rxn.name, assembly_id]
        else:
            x = ["reaction", rxn.name]
        return tuple(x + [tuple(inputs)] + [tuple(outputs)])

    # TODO: standard way to display SeqRecord

    # Create a DAG of the molecules and reactions
    g = nx.DiGraph()
    for q_i, (qk, result) in enumerate(design.results.items()):
        for a_i, a in enumerate(result.assemblies):
            dk = (qk, a_i)
            for r in a.reactions:
                r1 = rxn_key(r, a_i)
                g.add_node(r1, reaction=r, assembly_index=a_i, assembly=a)
                g.nodes[r1].setdefault("design_keys", set()).add(dk)
                for in_mol in r.inputs:
                    n1 = mol_key(in_mol)
                    g.add_node(n1, molecule=in_mol)
                    g.nodes[n1].setdefault("design_keys", set()).add(dk)
                    g.add_edge(n1, r1)
                for out_mol in r.outputs:
                    n2 = mol_key(out_mol)
                    g.add_node(n2, molecule=out_mol)
                    g.nodes[n2].setdefault("design_keys", set()).add(dk)
                    g.add_edge(r1, n2)
                    #
                    # if r.name == Reaction.Types.Assembly:
                    #     query_key_to_molecule_node

    # Validate the graph
    for n2, ndata2 in g.nodes(data=True):
        # validate each non-plasmid is expected to have at most one reaction
        molecule = ndata2.get("molecule", None)
        if molecule and molecule.type is not MoleculeType.types[C.PLASMID]:
            predecessors = list(g.predecessors(n2))
            if len(predecessors) > 1:
                msg = str(n2)
                for n1 in g.predecessors(n2):
                    msg += "\n\t" + str(n1)
                raise ValueError(
                    "Molecule {} has more than one reaction.\n{}".format(n2, msg)
                )

        reaction = ndata2.get("reaction", None)
        if reaction:
            predecessors = list(g.predecessors(n2))
            if not len(predecessors) == len(reaction.inputs):
                raise ValueError(
                    "The number of predecessors ({}) does not"
                    " match the number of inputs for the reaction"
                    " ({})".format(len(predecessors), len(reaction.outputs))
                )

            successors = list(g.successors(n2))
            if not len(successors) == len(reaction.outputs):
                raise ValueError(
                    "The number of successors ({}) does not"
                    " match the number of outputs for the reaction"
                    " ({})".format(len(successors), len(reaction.outputs))
                )

    # index the reactions and molecules in the graph
    # add the 'index' attribute to the node data
    reaction_node_dict = {}
    molecule_node_dict = {}
    sorted_nodes = list(nx.topological_sort(g))[::-1]
    for n in sorted_nodes:
        ndata = g.nodes[n]
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

    def used_in_designs(ndata):
        x = sorted(list(ndata["design_keys"]))
        return [{"design_key": _x[0], "assembly": _x[1]} for _x in x]

    property_reaction = []
    for n, index in reaction_node_dict.items():
        ndata = g.nodes[n]
        output_keys = sorted([mol_key(m) for m in ndata["reaction"].outputs])
        input_keys = sorted([mol_key(m) for m in ndata["reaction"].inputs])
        property_reaction.append(
            {
                "name": ndata["reaction"].name,
                "index": index,
                "used_in_assemblies": used_in_designs(ndata),
                "inputs": [molecule_node_dict[n] for n in input_keys],
                "outputs": [molecule_node_dict[n] for n in output_keys],
            }
        )

    # TODO: include information for SeqRecord?
    # TODO: convert SeqRecord to a JSON?
    property_molecule = []
    for n, index in molecule_node_dict.items():
        ndata = g.nodes[n]
        mol = ndata["molecule"]
        used_as_inputs = [reaction_node_dict[n] for n in g.successors(n)]
        used_as_outputs = [reaction_node_dict[n] for n in g.predecessors(n)]
        property_molecule.append(
            {
                "name": mol.type.name,
                "index": index,
                "sequence": str(mol.sequence.seq).upper(),
                "used_in_assemblies": used_in_designs(ndata),
                "used_as_input_to_reactions": used_as_inputs,
                "used_as_output_to_reactions": used_as_outputs,
            }
        )

    # TODO: in addition to status, add the following information about design/assemblies
    #       1. information about molecule and how it relates to the design
    #       2. how many fragments are used
    #       3. information about input parameters
    output_json = {
        "designs": design.status,
        "reactions": property_reaction,
        "molecules": property_molecule,
    }

    print(json.dumps(output_json, indent=2))
    # for n2 in g.nodes():
    #     print(n2)
    #     for n1 in g.predecessors(n2):
    #         print('\t' + str(n1))
    # e1 = list(g.edges())
    # e2 = set(g.edges())
    # print(e1)
    # print(len(e1))
    # print(len(e2))

    # for n1, n2, edata in g.edges(data=True):
    #     print(n1, n2)
