# DASi DNA Design

**DASi** is an automatic DNA cloning plan designer aimed for operating on small budgets
by focusing on material re-use.

The software converts a nucleotide sequence, or a library of sequences, to an executable
 molecular assembly plan while
optimizing material cost, assembly efficiency, and assembly time.

The software goals are reminiscent of j5 or Teselegen but focused on:
1. having a dead-simple user interface and
1. utilizing information about current laboratory inventory in its optimization
algorithm.

### Planned Features

* Golden-gate support
* heirarchical assembly
* library support (with bayesian search to optimize shared parts)
* front-end
* connection to fabrication facility

### Use cases

* developing cloning plans from computer-generated sequences
* developing cloning plans for human-generated sequences
* developing plans for users that do not know the intricacies of molecular biology

### Other related repos used in this project:

* pyblastbio - python BLAST wrapper
* primer3-py-plus - python wrapper around Primer3
* loggable-jdv - logging class
* benchlingapi - Python BenchlingAPI
