from shoestring import AlignmentContainer, Constants


def test_load(new_bio_blast, new_primer_blast):
    blast = new_bio_blast()
    blast.quick_blastn()

    primer_blast = new_primer_blast()
    primer_blast.quick_blastn_short()

    results = blast.get_perfect()
    primer_results = primer_blast.get_perfect()

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

    products = container.groups_by_type[Constants.PCR_PRODUCT]
