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

**Possible Failure Modes**

* overhang are homologous to other overhangs. This could be translated to some global score.
* Efficiency of PCR amplification could be low. This could be determined by Primer3.
* PCR amplification could have multiple products. This could be determined by Primer3.

**Multi-plasmid optimization**

0. Find shared edges between assembly graphs, we concatenate all of the queries separated by ???'s and align the subjects.
We then break down subject alignments to each original query.
1. We also need to find shared primers somehow. Enumerating every shared primer would be very comp intensive in the
same way we do with shared edges is probably not possible.
2. After we align primers to each query

1. For each plasmid, enumerate all possible shared edges. For example, if a assembly graph has a single shared
edge with another assembly graph, then there are two solutions, one in which the shared edge is not used and 
the other in which the shared edge is used (A->B cost /= 2). In another example, G1 shares edges with G2 and G3.
Then there are 3 solutions, not sharing an edge (A-> cost /= 1). Sharing a single edge (A->B cost /= 2) and sharing
both edges (A->B cost /= 3). The total number of solutions will be equal to `s1*s2*s3*...` where `si` is the maximum
number of shared edges for edge `i`. We should optionally provide a way to boost solutions that share edges.
2. Once all solutions have been enumerated with corresponding costs, we find compatible solutions (those with shared edges)
and simply choose the smallest cost.

**Hierarchical optimization**

How do you optimize a plasmid AND then potentially use that as a template for the next plasmid?

## Bug Fixes

1. `check_design` method in testing needs to be updated, since its returning np.inf.
2. plot assembly using matplotlib

**Primer Optimization**


* Find primers that have high sequence similarity. Make internal edge cost high in those cases.
* Primers with different Tms. Internal edge efficiency.
* Assembly junctions with high sequence similarity.