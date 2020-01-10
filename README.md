# DASi DNA Design

[![PyPI version](https://badge.fury.io/py/dasi.svg)](https://badge.fury.io/py/dasi)

**DASi** is an automatic DNA cloning plan designer aimed for operating on small budgets
by focusing on material re-use.

The software converts a nucleotide sequence, or a library of sequences, to an executable
 molecular assembly plan while optimizing material cost, assembly efficiency, and assembly time.

The key design paradigm for DASi is that *no molecular biology expertise is required to use DASi*. Complete novices should be able to use the software to design and construct new genetic sequences.

The software goals are reminiscent of j5 or Teselegen but focused on:
1. having a dead-simple user interface and
1. utilizing information about current laboratory inventory in its optimization
algorithm.

### Status

DASi is currently under development.

### Usage

DASi completely automates the cloning design work, finding approximately optimal solutions for cloning steps, preferentially using existing plasmids, linear DNA fragments, and primers to design semi-optimal cloning steps and designs.

The following command designs the cloning steps for a library of designs. The user only needs to specify the sequences they wish to construct and currently available primers, DNA templates, and linear DNA as *.genbank* or *.fasta* files. DASi handles all design aspects. *No molecular biology expertise is required to use DASi.*

```bash
dasi library_design --designs mydesigns/*.gb --fragments fragments/*.gb --primers primers.fasta --templates plasmids/*.gb --cost_model cost.b --out results
```

#### Customization

DASi optimization parameters are completely customizable. The following

* primer synthesis costs
* primer design parameters
* synthetic fragment costs
* vendor-specific synthetic fragment complexity
* sequence dependent plasmid assembly efficiencies
* optimizing over efficiency vs material costs

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
