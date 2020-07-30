===============
Getting Started
===============

**Python Installation**

.. code-block:: bash

    pip install dasi
    pyblast install

**Using DASi as a Docker command**

.. code-block:: bash

    # TODO: put the Docker command here

**Hello World**

.. code-block:: python

    import dasi
    import json

    # create a fake database and create 10 fake designs
    design = dasi.Design.fake(n_design=10)
    design.run()
    print(json.dumps(design.out(), indent=2))

**Getting Help**

.. code-block:: python

    dasi.Design.help()
    dasi.Design.documentation()



