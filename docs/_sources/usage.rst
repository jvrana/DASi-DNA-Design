Usage
=====

Getting started
---------------

Installation
^^^^^^^^^^^^

Install or upgrade using pip3

.. code-block:: bash

    pip install dasi -U

Running examples
^^^^^^^^^^^^^^^^

You can produce and run a randomized design by using the following:

.. code-block:: python

    from dasi import Design
    import json

    design = Design.fake(n_designs=1)
    design.run(n_paths=1)
    print(json.dumps(design.out(), indent=2))

Advanced
--------

Using inventory information with DASi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Additional sequence information can be added to sequences before DASi design,
such that this information appears in the final result. For example, one
may wish to add inventory information to specific primers or templates
so you can find it in the lab after DASi finishes its design.

Lets start with some sequences and an empty Design

.. code-block:: python

    from dasi import Design

    design = Design()

We open our sequences as BioPythons SeqRecord objects.

.. code-block:: python

    from Bio import SeqIO
    from glob import glob

    primers = []
    with primer_path in glob("primers/*.gb"):
        primers.append(SeqIO.read(primer_path, format='genbank'))

    fragments = []
    with primer_path in glob("fragments/*.gb"):
        primers.append(SeqIO.read(primer_path, format='genbank'))

    # etc.

We add additional information to the annotations.

.. code-block:: python

    for f in fragments:
        f.annotations['location'] = 'benchtop'


The annotations should appear in the results
(e.g. `results['molecules'][0]['customFields']`)

.. code-block:: python

    design.add_material(fragments=fragments, plasmids=plasmids, primers=primers, queries=queries)
    design.run()
    print(design.out())


Adjusting design parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Documentation coming soon.
