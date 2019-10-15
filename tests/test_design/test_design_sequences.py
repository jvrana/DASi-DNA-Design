import json
import os
from typing import Dict
from typing import List

import dill
import jdna
import pytest
from primer3plus.utils import anneal
from primer3plus.utils import reverse_complement as rc
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.design import Design
from dasi.models import Reaction
from dasi.models.assembly import design_edge
from dasi.models.assembly import design_primers
from dasi.utils import Region

gfp = "ATGGTCTCTAAGGGTGAAGAATTGTTCACCGGTGTCGTCCCAATCTTGGTCGAATTGGACGGGGACGTCAACGGTCACAAGTTCTCTGTCTCTGGTGAAGGTGAAGGTGACGCTACCTACGGTAAGTTGACCTTGAAGTTCATCTGTACCACCGGTAAGTTGCCAGTCCCATGGCCAACCTTGGTCACCACCTTCGGTTACGGTGTCCAATGTTTCGCTAGATACCCAGACCACATGAAGCAACACGACTTCTTCAAGTCTGCTATGCCAGAAGGTTACGTCCAAGAAAGAACCATCTTCTTCAAGGACGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTCGAAGGTGACACCTTGGTCAACAGAATCGAATTGAAGGGTATCGACTTCAAGGAAGACGGTAACATCTTGGGTCACAAGTTGGAATACAACTACAACTCTCACAACGTCTACATCATGGCTGACAAGCAAAAGAACGGTATCAAGGTCAACTTCAAGATCAGACACAACATCGAAGACGGTTCTGTCCAATTGGCTGACCACTACCAACAAAACACCCCAATCGGTGACGGTCCAGTCTTGTTGCCAGACAACCACTACTTGTCTACCCAATCTGCTTTGTCTAAGGACCCAAACGAAAAGAGAGACCACATGGTCTTGTTGGAATTCGTCACCGCTGCTGGTATCACCCACGGTATGGACGAATTGTACAAGTAA"


def test_region_invert():
    """Check Region.invert.

    This is essential to the primer design algorithm.
    """
    r = Region(10, 20, 50, cyclic=True)
    r1, r2 = r.invert()
    assert r2 is None
    assert list(r1) == list(range(20, 50)) + list(range(10))


def test_primer_design():
    region = Region(100, 300, len(gfp), cyclic=True)
    pairs, explain = design_primers(gfp, region, None, None)
    print(json.dumps(pairs, indent=1))

    for pair in pairs.values():
        assert pair["LEFT"]["location"][0] == 100
        assert pair["RIGHT"]["location"][0] == 299


def test_primer_design2():
    i = len(gfp) - 50
    j = 100
    region = Region(i, j, len(gfp), cyclic=True)
    pairs, explain = design_primers(gfp, region, None, None)
    # print(json.dumps(pairs, indent=1))
    for pair in pairs.values():
        print(pair["LEFT"]["location"])
        assert pair["LEFT"]["location"][0] == i
        assert pair["RIGHT"]["location"][0] == j - 1


def test_primer_design_overorigin():

    template = "AGCTGGAGAATTGCCATGTAGATGTTCATACAATCGTCAAATCATGAAGGCTGGAAAAGCCCTCCAAGATCCCCAAGACCAACCCCAACCCACCCACCGTGCCCACTGGCCATGTCCCTCAGTGCCACATCCCCACAGTTCTTCATCACCTCCAGGGACGGTGACCCCCCCACCTCCGTGGGCAGCTGTGCCACTGCAGCACCGCTCTTTGGAGAAGGTAAATCTTGCTAAATCCAGCCCGACCCTCCCCTGGCACAACGTAAGGCCATTATCTCTCATCCAACTCCAGGACGGAGTCAGTGAGGATGGGGCTCTAGCTCTAGAGCTTGATATCGAATTCCTGCAGCCCCGGGACAGCCCCCCCCCAAAGCCCCCAGGGATGTAATTACGTCCCTCCCCCGCTAGGGGGCAGCAGCGAGCCGCCCGGGGCTCCGCTCCGGTCCGGCGCTCCCCCCGCATCCCCGAGCCGGCAGCGTGCGGGGACAGCCCGGGCACGGGGAAGGTGGCACGGGATCGCTTTCCTCTGAACGCTTCTCGCTGCTCTTTGAGCCTGCAGACACCTGGGGGGATACGGGGAAAAAGCTTTAGGCTGAAAGAGAGATTTAGAATGACAGAATCATAGAACGGCCTGGGTTGCAAAGGAGCACAGTGCTCATCCAGATCCAACCCCCTGCTATGTGCAGGGTCATCAACCAGCAGCCCAGGCTGCCCAGAGCCACATCCAGCCTGGCCTTGAATGCCTGCAGGGATGGGGCATCCACAGCCTCCTTGGGCAACCTGTTCAGTGCGTCACCACCCTCTGGGGGAAAAACTGCCTCCTCATATCCAACCCAAACCTCCCCTGTCTCAGTGTAAAGCCATTCCCCCTTGTCCTATCAAGGGGGAGTTTGCTGTGACATTGTTGGTCTGGGGTGACACATGTTTGCCAATTCAGTGCATCACGGAGAGGCAGATCTTGGGGATAAGGAAGTGCAGGACAGCATGGACGTGGGACATGCAGGTGTTGAGGGCTCTGGGACACTCTCCAAGTCACAGCGTTCAGAACAGCCTTAAGGATAAGAAGATAGGATAGAAGGACAAAGAGCAAGTTAAAACCCAGCATGGAGAGGAGCACAAAAAGGCCACAGACACTGCTGGTCCCTGTGTCTGAGCCTGCATGTTTGATGGTGTCTGGATGCAAGCAGAAGGGGTGGAAGAGCTTGCCTGGAGAGATACAGCTGGGTCAGTAGGACTGGGACAGGCAGCTGGAGAATTGCCATGTAGATGTTCATACAATCGTCAAATCATGAAGGCTGGAAAAGCCCTCCAAGATCCCCAAGACCAACCCCAACCCACCCACCGTGCCCACTGGCCATGTCCCTCAGTGCCACATCCCCACAGTTCTTCATCACCTCCAGGGACGGTGACCCCCCCACCTCCGTGGGCAGCTGTGCCACTGCAGCACCGCTCTTTGGAGAAGGTAAATCTTGCTAAATCCAGCCCGACCCTCCCCTGGCACAACGTAAGGCCATTATCTCTCATCCAACTCCAGGACGGAGTCAGTGAGGATGGGGCTCTAGCGGGGGGATCCGATGTCGACACGCGTGCATGCGCCGATACGAAGGTTTTCTCCAGCGAAGGTCGGGCAGGAAGAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGGGAACGTGATTGAATAACTTTGGCCTCGACTCTGTCAACTGACTTCCCCCGTCGTTCACTGCCGTATAGGCAGCATCTTTAGAATAGCTCAGAGGCCGAGGGTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCATTCCTGTTCACTGCCGTATAGGCAGCCCTTTATCTCCCACGTGCGCTTTCTCCCTTCTCCTTTTTTCTAGGCTTCAATAAAGGAGCGAGCACCCGTGCCGGAGACCCACAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGATCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGGAAGACCCAATGGTCGGCGGGACCAGGGAGTTTAAACTAGCATCGCGATAAGCTCTAGAGGGACAGCCCCCCCCCAAAGCCCCCAGGGATGTAATTACGTCCCTCCCCCGCTAGGGGGCAGCAGCGAGCCGCCCGGGGCTCCGCTCCGGTCCGGCGCTCCCCCCGCATCCCCGAGCCGGCAGCGTGCGGGGACAGCCCGGGCACGGGGAAGGTGGCACGGGATCGCTTTCCTCTGAACGCTTCTCGCTGCTCTTTGAGCCTGCAGACACCTGGGGGGATACGGGGAAAAAGCTTTAGGCTGAAAGAGAGATTTAGAATGACAGAATCATAGAACGGCCTGGGTTGCAAAGGAGCACAGTGCTCATCCAGATCCAACCCCCTGCTATGTGCAGGGTCATCAACCAGCAGCCCAGGCTGCCCAGAGCCACATCCAGCCTGGCCTTGAATGCCTGCAGGGATGGGGCATCCACAGCCTCCTTGGGCAACCTGTTCAGTGCGTCACCACCCTCTGGGGGAAAAACTGCCTCCTCATATCCAACCCAAACCTCCCCTGTCTCAGTGTAAAGCCATTCCCCCTTGTCCTATCAAGGGGGAGTTTGCTGTGACATTGTTGGTCTGGGGTGACACATGTTTGCCAATTCAGTGCATCACGGAGAGGCAGATCTTGGGGATAAGGAAGTGCAGGACAGCATGGACGTGGGACATGCAGGTGTTGAGGGCTCTGGGACACTCTCCAAGTCACAGCGTTCAGAACAGCCTTAAGGATAAGAAGATAGGATAGAAGGACAAAGAGCAAGTTAAAACCCAGCATGGAGAGGAGCACAAAAAGGCCACAGACACTGCTGGTCCCTGTGTCTGAGCCTGCATGTTTGATGGTGTCTGGATGCAAGCAGAAGGGGTGGAAGAGCTTGCCTGGAGAGATACAGCTGGGTCAGTAGGACTGGGACAGGC"
    region = Region(2020, 1616, len(template), cyclic=True)

    adjusted_template = region.get_slice(template) + region.invert()[0].get_slice(
        template
    )

    rprimer = "CGCTGGAGAAAACCTTCGTATCGGCgcatgcacgcgtgtcgacatcg"

    assert (
        rc("CGCTGGAGAAAACCTTCGTATCGGCgcatgcacgcgtgtcgacatcg".upper())
        in adjusted_template.upper()
    )

    fwd, rev = anneal(adjusted_template, [rprimer])
    print(rev[0]["top_strand_slice"])
    # reindexed = region.get_slice(template) + region.invert()[0].get_slice(
    #     template
    # )
    # assert expected == reindexed

    rprimer = "CGCTGGAGAAAACCTTCGTATCGGCgcatgcacgcgtgtcgacatcg"
    pairs, explain = design_primers(template, region, None, rseq=rprimer)
    print(json.dumps(pairs, indent=1))
    print(explain)
    assert pairs


class TestExpectedSequences:
    def design_for_assembly(self, design, assembly):
        reactions = []
        for n1, n2, edata in assembly.edges():
            reaction = design_edge(assembly, n1, n2, design.seqdb)
            if reaction:
                reactions.append(reaction)
        return reactions

    def validate_pcr(self, reaction: Reaction, length_only=False):
        primers = [m for m in reaction.inputs if m.type.name == "PRIMER"]
        template = [m for m in reaction.inputs if m.type.name == "TEMPLATE"][0]
        template_seq = jdna.Sequence(
            template.sequence, cyclic=template.alignment_group.subject_region.cyclic
        )
        if template.alignment_group.subject_region.direction == -1:
            template_seq = template_seq.reverse_complement()
        primer_seqs = [jdna.Sequence(p.sequence, cyclic=False) for p in primers]

        product = jdna.Reaction.pcr(template_seq, primer_seqs)[0]

        if length_only:
            assert len(product) == len(reaction.outputs[0].query_region)
            return
        assert str(product).upper() == str(reaction.outputs[0].sequence.seq).upper()

    def reactions_to_assembly(self, reactions):
        seqs = []
        for r in reactions:
            sequence = r.outputs[0].sequence
            seqs.append(jdna.Sequence(sequence, cyclic=False))
        assemblies = jdna.Reaction.cyclic_assemblies(seqs)
        return assemblies

    # @pytest.fixture(scope="module")
    # def processed_assemblies(self, reactions_dict):
    #     d = {}
    #     for qk, rlist in reactions_dict.items():
    #         d[qk] = []
    #         for reactions in rlist:
    #             d[qk].append(self.reactions_to_assembly(reactions))
    #     return d

    @pytest.fixture(scope="module")
    def reactions_dict(
        self, multi_processed_results
    ) -> Dict[str, List[List[Reaction]]]:
        design, results = multi_processed_results
        reaction_dict = {}
        for qk, result in results.items():
            reaction_dict[qk] = []
            for assembly in result.assemblies:
                reaction_dict[qk].append(self.design_for_assembly(design, assembly))
        return reaction_dict

    def test_check_pcr_product(self, multi_processed_results, reactions_dict):
        """Tests expected product length and/or product sequence."""
        tested = False
        for qk, reactions_list_of_lists in reactions_dict.items():
            for reactions_list in reactions_list_of_lists:
                for r in reactions_list:
                    if r.name == "PCR":
                        self.validate_pcr(r)
                        tested = True
        assert tested

    def test_num_assemblies(self, reactions_dict):
        for qk, rlist in reactions_dict.items():
            for reactions in rlist:
                assemblies = self.reactions_to_assembly(reactions)
                assert len(assemblies) == 2

    def test_call_reactions_from_assembly(self, multi_processed_results):
        design, results = multi_processed_results
        for result in results.values():
            for assembly in result.assemblies:
                assert assembly.reactions

                # print("REACTIONS")
                # for r in reactions:
                #     print(r)
        #         assemblies = self.reactions_to_assembly(reactions)
        #         print("ASSEMBLIES")
        #         for a in assemblies:
        #             print(a)

    # def test_has_gibson_assembly(self, processed_assemblies):
    #     failed = []
    #     for qk, assemblies in processed_assemblies.items():
    #         if not assemblies:
    #             failed.append(qk)
    #     assert not failed
    #
    # def test_has_only_one_assembly(self, processed_assemblies):
    #     failed = []
    #     i = 0
    #     for qk, list_of_assemblies in processed_assemblies.items():
    #         for ia, assemblies in enumerate(list_of_assemblies):
    #             if len(assemblies) != 2:
    #                 print(assemblies)
    #                 failed.append((i, ia, qk, len(assemblies)))
    #         i += 1
    #     assert not failed
    #
    # def test_gibson_assembly_product_size(self, processed_assemblies,
    #                                       multi_processed_results):
    #     design, _ = multi_processed_results
    #     for qk, assemblies in processed_assemblies.items():
    #         expected_record = design.seqdb[qk]
    #         product = assemblies[0].product
    #         assert len(product) == len(expected_record)

    # @pytest.mark.parametrize('stringency', [0, 1, 2])
    # def test_gibson_assembly(self, multi_processed_results, reactions_dict, stringency):
    #     design, results = multi_processed_results
    #     for qk, reaction_list_of_lists in reactions_dict.items():
    #         assemblies = results[qk].assemblies
    #         for i, reaction_list in enumerate(reaction_list_of_lists):
    #             assembly = assemblies[i]
    #             print(assembly.to_df())
    #             print("*" * 100)
    #             print("REACTIONS:")
    #             for r in reaction_list:
    #                 print(r)
    #
    #             check_num_assemblies = False
    #             if stringency == 0:
    #                 check_num_assemblies = True
    #                 self.validate_assembly(reaction_list, check_num_assemblies=True)
    #             elif stringency ==
