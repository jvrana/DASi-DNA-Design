**10/10/19**

*interface*
* add a basic CLI
* expose parameters to input file
* convert assembly to jdna and print assembly to print to file.

*logging/debugging*
* add a level inbetween INFO and ERROR for logger for status update to CLI
* add log to text file

*Sequence Database*
* save filename and path for SeqRecord on load

*Coalition optimization*
* shared templates? select from same templates, if possible
* shared sequences and fragments?
* shared primers?
* AlignmentGroups should be shared directly. When a sequence is designed,
this same instance should be shared among other members of the coalition.

*Sequences*
* [IN PROGRESS] costs should be saved to Reactions
* self blast to eliminate potential overhangs for mis-assembly
* max size lim and tests?

*Tests*
* add tests for checking gibson assembly and PCR products
* [IN PROGRESS] tests for size limitations

*Refactoring*
* rename `sequence_design` to `sequence_designer`
* rename `design_algorithms` to `design_optimization` or `utils`
* replace jdna with something faster

**Other**
* Resolve TODOs
* Basic design of sequences
* Bayesian optimization to design junctions
* Back convert assemblies to a plan
* better cost estimates
* break apart long contigs
* primer3 design
* after finding shortest paths, refine weights using simulation, update weights, repeat until there is no improvements
* force parts
* avoid repeats
* consider design flexibility in optimization algorithm
* primers for amplification of synthesis fragments
* contrived tests with primer and gene designs

*algorithm improvements*
* eliminating non-sensical cycles. (e.g. 1000 -> 500 -> 501 -> 1000)
* black-listed homology edges during modified floyd-warshall algorithm.
During the computation of shortest path lengths between pairs. Blacklisted
edges maintained separately. Would need keep an up-to-date matrix
of blacklisted node as in floyd-warshall. Considered path would
be normal cost + blacklisted (0 or inf). Blacklisted edges determined
by blast searching goal plasmids with themselves and finding alignment
groups whose ends fall within these blacklisted region.


**Possible Failure Modes**

* overhang are homologous to other overhangs. This could be translated to some global score.
* Efficiency of PCR amplification could be low. This could be determined by Primer3.
* PCR amplification could have multiple products. This could be determined by Primer3.

**Multi-plasmid optimization**

1. Perform design-on-design blast search. GET <BLAST_RESULTS_JSON>
2. Remove self alignments from blast <BLAST_RESULTS_JSON>
3. Group alignments in a un-directed graph

**Hierarchical optimization**

How do you optimize a plasmid AND then potentially use that as a template for the next plasmid?

## Bug Fixes

1. `check_design` method in testing needs to be updated, since its returning np.inf.
2. plot assembly using matplotlib

**Primer Optimization**


* Find primers that have high sequence similarity. Make internal edge cost high in those cases.
* Primers with different Tms. Internal edge efficiency.
* Assembly junctions with high sequence similarity.

**Gap partitioning**

Partition Score:
https://colab.research.google.com/drive/19YZ7Jnx-4cfJAbNoDF7wgMheCqRCehPZ