Tests
=====

Tests are written to run with pytest.

To run tests, first install trident's dependencies:

.. code::

    make

To run tests, run:

.. code::

    make tests

To run doctests located in `user/examples`.

.. code::

    make doctest

Test coverage
-------------

Covering all of the models and endpoints for Trident is very difficult.
Tests coverage is not 100%.

For the most part, the `Marshaller`, `utilities`, `aqhttp`, `aqsession`, and
`baseclass` have near 100% coverage.

For testing of specific Aquarium models, tests are found in
'tests/test\_models/models.' Many of these model tests use 'monkeypatch'
to *intercept* requests and return expected information.
Writing these tests take a long time and so not all model tests are comprehensive.

Request recording
~~~~~~~~~~~~~~~~~

By default, live requests are recorded automatically via
`VCRpy <https://vcrpy.readthedocs.io/en/latest/installation.html>`_

Information about each request is stored in a **fixtures/vcr_cassettes/*.yaml**
and retrieved on a as-needed basis.


Testing with live and fake connections
--------------------------------------

In 'tests/conftest.py' two fixtures are defined `fake_session` and
`session.` `fake_session` contains a fake AqSession in which the
login information has been faked. `session` contains a potentially
live connection to Aquarium using login information located in
'tests/secrets/config.json.secret'

Though we should not run tests against a production or nursery server,
to run live tests with an Aquarium server, replace the
`tests/secrets/config.json.secret.template` and a new
`tests/secrets/config.json.secret` containing your login information

.. code-block:: JSON

    {
        "login": "yourlogin",
        "password": "password",
        "aquarium_url": "http://aquarium.org"
    }
