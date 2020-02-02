import json

import networkx as nx
import pylab as plt

from dasi import Design
from dasi import LibraryDesign


def test_output():
    design = Design.fake(n_designs=3, n_primers_from_templates=1000)
    design.run()

    def mol_key(mol1, assembly_index):
        return (mol1.type.name, assembly_index, str(mol1.sequence.seq).upper())

    def add_edge(g, a_i, mol1, mol2):
        n1 = mol_key(mol1, a_i)
        n2 = mol_key(mol2, a_i)
        g.add_node(n1, molecule=mol1)
        g.add_node(n2, molecule=mol2)
        g.add_edge(n1, n2)

    g = nx.DiGraph()

    molecules = {}
    reactions = []
    assemblies = []
    for q_i, (qk, result) in enumerate(design.results.items()):
        for a_i, a in enumerate(result.assemblies):
            assemblies.append({"query_key": qk, "reactions": []})
            for r in a.reactions:
                for mol in r.inputs + r.outputs:
                    k = mol_key(mol, a_i)
                    if k not in molecules:
                        molecules[k] = {
                            "index": None,
                            "type": mol.type.name,
                            "sequence": str(mol.sequence.seq).upper(),
                        }
                        molecules[k]["index"] = len(molecules)
                reactions.append(
                    {
                        "index": len(reactions),
                        "name": r.name,
                        "assembly": a_i,
                        "query": qk,
                        "inputs": [
                            molecules[mol_key(m, a_i)]["index"] for m in r.inputs
                        ],
                        "outputs": [
                            molecules[mol_key(m, a_i)]["index"] for m in r.outputs
                        ],
                    }
                )
                assemblies[-1]["reactions"].append(reactions[-1]["index"])
                #     # for out_mol in r.outputs:
                #     #     add_edge(g, q_i, a_i, in_mol, out_mol)

    molecules = [ndata["molecule"] for _, ndata in g.nodes(data=True)]

    output_json = {
        "assemblies": assemblies,
        "molecules": sorted(molecules, key=lambda x: x["index"]),
        "reactions": reactions,
    }

    # for r in reactions:
    #     output_json['reactions'].append(
    #         {
    #             'name': r.name
    #         }
    #     )

    print(json.dumps(output_json, indent=2))

    # nx.draw(g)
    # plt.show()
