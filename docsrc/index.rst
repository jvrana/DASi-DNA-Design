Introduction
------------

DASi is an automatic DNA cloning plan designer aimed for operating on small budgets
by focusing on material re-use.

The software converts a nucleotide sequence, or a library of sequences,
to an executable molecular assembly plan while optimizing material cost,
assembly efficiency, and assembly time.


User Documentation
------------------

The user documentation contains high-level information for users of DASi.

.. toctree::
   :maxdepth: 1

   usage
   schemas/schemas.rst

API Reference
-------------

Detailed information about classes and methods of DASi.

.. toctree::
   :maxdepth: 1

   design
   models
   cost
   command_line
   utils
   exceptions
   schemas
   constants

Developer Documentation
-----------------------

The developer documentation conatins information on how to contribute to DASi.

.. toctree::
   :maxdepth: 1

   guidelines

Features and Improvements Road Map
----------------------------------

* ^0.0.1
    - Cleaner JSON output
    - Scoring PCR reactions for primer misprimings
    - Scoring assembly efficiencies
    - Logging to file
    - PDF design reports

* ^0.1
    - faster execution times
    - improved design inputs
        + selection of synthesis vendors
        + primer + PCR cost
        + improved API for adding materials
        + control over design algorithms and parameters

* ^0.2
    - linear fragment design
        + support for assembly by primer annealing

* ^0.3
    - improved outputs
        + SBOL (XML) output
        + graph output

* ^1.0
    - deployment
        + web interface
        + dockerized container

* ^2.0
    - Additional assembly strategies
        + Golden-gate support
        + Linear assemblies
        + Assembly by primers

* ^3.0
    - Hierarchical clonal assembly

