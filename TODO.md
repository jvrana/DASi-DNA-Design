# TODO (11/21/19)

*Sequence Database*
* save filename and path for SeqRecord on load

*better way to add notes to graph edge*

*graph post-processing config*

## Coalition optimization
* Make a new class called DesignConfig
    + Should contain min_primer_anneal. This value should be used
        for calculating primer3 design parameters (as in primer3plus.utils.anneal).
        This same value minus one, should be used to find misprimings in the Graph
        post processor. User should be able to relax or tighten the post processor.
    + Could contain n_jobs, n_paths
    + Other post processor parameters
    + cost parameters?
* shared templates? select from same templates, if possible
* shared sequences and fragments?
* shared primers?

## IDT Synthesis Complexity

### Implementing Complexity Scores for Gaps

#### Milestones

0. Implement ability to partition B->A edges. DONE.
1. Implement complexity to efficiency score for all gapped plasmids

**Other**
* Resolve TODOs
* Basic design of sequences
* avoiding homology repeats
* sharing primer designs across multiple designs
* sharing fragments across multiple designs
* add a level inbetween INFO and ERROR for logger for status update to CLI
* add log to text file

**Gap partitioning**

Partition Score:
https://colab.research.google.com/drive/19YZ7Jnx-4cfJAbNoDF7wgMheCqRCehPZ

**