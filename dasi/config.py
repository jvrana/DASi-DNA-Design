class Config:
    """Global DASi design configuration."""

    PRIMER3_MIN_ANNEAL_CHECK = 12  #: min n bases to check for misprimings
    PRIMER3_N_RUNS = (
        15  #: number of optimization steps to use for primer3 primer designs
    )
    ASSEMBLY_COST_THRESHOLD = 10000  #: maximum cost of an assembly

    ############
    # PARAMETERS
    ############
    PRIMER_MIN_BIND = 14  #: minimum primer binding for searching for primer alignments
    MIN_OVERLAP = 15  #: minimum overlap for searching for overlapping alignments
    MAX_HOMOLOGY = 100  #: maximum overlap for searching for overlapping alignments
    INF = 10.0 ** 6  #: almost infinity

    class SequenceScoringConfig:
        """Configuration for scoring sequences."""

        stats_repeat_window = 14  #: length of kmer to find sequence repeats
        stats_window = 20  #: length of window for sliding window calculations
        stats_hairpin_window = 20  #: length of kmer to find sequence hairpins
        mispriming_min_anneal = 12  #: minimum bp to look for misprimings
        mispriming_max_anneal = 30  #: maximum expected bp for a primer
        mispriming_penalty = 0.5  #: multiplier to apply to each mispriming

        #: threshold for sequence complexity that is not synthesizable by a vendor
        complexity_threshold = 10

        #: edges with efficiencies below this threshold will be removed.
        edge_threshold = 0.05

        #: efficiency value for sequence that is not synthesizable
        not_synthesizable_efficiency = 0.1

        #: multiplier to apply to the efficiency if PCR is within the given length range
        pcr_length_range_efficiency_multiplier = [
            (4000, 5000, 0.8),
            (5000, 6000, 0.5),
            (6000, 1000, 0.2),
        ]

        partition_overlap = 30  #: overlap when partitioning a highly complex sequence
        partition_step_size = 10  #: step size to approximate optimal partition
        SCORE_COMPLEXITY = "SCORE_COMPLEXITY"
        SCORE_MISPRIMINGS = "SCORE_MISPRIMINGS"
        SCORE_LONG = "SCORE_LONG_PCR_PRODUCTS"
        PARTITION = "PARTITION_SYNTHETIC_SEQUENCES"
        post_process_stages = (
            SCORE_COMPLEXITY,
            SCORE_LONG,
            SCORE_MISPRIMINGS,
            PARTITION,
        )
