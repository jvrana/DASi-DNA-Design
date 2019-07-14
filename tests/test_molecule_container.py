from shoestring import AlignmentContainer, Constants


def test_load(new_bio_blast, new_primer_blast):
    blast = new_bio_blast()
    blast.quick_blastn()

    primer_blast = new_primer_blast()
    primer_blast.quick_blastn_short()

    results = blast.get_perfect()
    primer_results = primer_blast.get_perfect()

    def perfect_subject(data):
        if data['strand'] == 1 and data['start'] == 1 and data['end'] == data['origin_sequence_length']:
            return True
        elif data['strand'] == -1 and data['end'] == 1 and data['start'] == data['origin_sequence_length']:
            return True

    primer_results = [p for p in primer_results if perfect_subject(p['subject'])]
    print("Number of perfect primers: {}".format(len(primer_results)))
    # primer_results = [p for p in primer_results if p['subject']['start'] == 1]

    # combine the sequence databases (for now)
    seqdb = {}
    seqdb.update(blast.seq_db.records)
    seqdb.update(primer_blast.seq_db.records)

    container = AlignmentContainer(seqdb)

    # load the results (potential PCR Products)
    container.load_blast_json(results, Constants.PCR_PRODUCT)

    # load the primer results
    container.load_blast_json(primer_results, Constants.PRIMER)

    # group by query_regions
    groups = container.alignment_groups

    print("Number of types: {}".format(len(container.groups_by_type)))
    print("Number of groups: {}".format(len(groups)))

    G = container.build_assembly_graph()
    assert G.number_of_edges()
