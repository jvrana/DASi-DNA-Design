# DASi DNA Design

[![PyPI version](https://badge.fury.io/py/dasi.svg)](https://badge.fury.io/py/dasi)

**DASi** is an automatic DNA cloning plan designer aimed for operating on small budgets
by focusing on material re-use.

The software converts a nucleotide sequence, or a library of sequences, to an executable
 molecular assembly plan while optimizing material cost, assembly efficiency, and assembly time.

The key design paradigm for DASi is that *no molecular biology expertise* is required to use DASi. Complete novices should be able to use the software to design and construct new genetic sequences. This also enables automated software programs to automatically design and construct new genetic sequences.

The software goals are reminiscent of j5 or Teselegen but focused on:
1. A dead-simple API usable by lab novices, experts or automated software programs.
1. Utilizing information about current laboratory inventory in its optimization
algorithm to minimize costs and turn-around time

### Status

DASi is currently under development funded by the DARPA Synergistic Discovery and Design program. DASi is currently being used to connect automatically generate DNA designs to automated biological fabrication facilities (e.g. University of Washington Biofab).

### Usage

DASi completely automates the cloning design work, finding approximately optimal solutions for cloning steps, preferentially using existing plasmids, linear DNA fragments, and primers to design semi-optimal cloning steps and designs.

The following command designs the cloning steps for a library of designs. The user only needs to specify the sequences they wish to construct and currently available primers and DNA templates as *.genbank* or *.fasta* files. DASi handles all design aspects. *No molecular biology expertise is required to use DASi.*

```bash
dasi library_design --designs mydesigns/*.gb --fragments fragments/*.gb --primers primers.fasta --templates plasmids/*.gb --cost_model cost.b --out results
```

#### Customization

DASi optimization parameters are completely customizable. The following are examples of parameters and aspects of DASi that are customizable:

* primer synthesis costs
* primer design parameters
* synthetic fragment costs
* vendor-specific synthetic fragment complexity
* sequence dependent plasmid assembly efficiencies
* optimizing over efficiency vs material costs
* etc.

### Planned Features

* Golden-gate support
* heirarchical assembly
* library support (with bayesian search to optimize shared parts)
* front-end
* connection to fabrication facility

### DASi optimization problem

Briefly, DASi approximates a solution the following optimization problem: 

```Given a set of 'goal' double-stranded sequences, a set of available single-stranded and double-strand sequences, and a set of actions that can create new sequences, find the optimal set of operations that produces the 'goal' sequences.```

Formalization of this optimization problem is coming soon.
