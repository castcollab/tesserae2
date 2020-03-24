Overview of Tesserae2
=====================

Tesserae2: Fast recombination-aware global and local alignment.

Current version: 1.0.0

Audience
--------


Free software
-------------

Tesserae2 is free software; you can redistribute it and/or modify it under the
terms of the Apache License version 2.0.  Contributions are welcome. Please join us on `GitHub <https://github.com/winni2k/cortexpy>`_.


Installation
------------

::

    pip install .


Documentation
-------------


Bugs
----

For bugs, please raise a `GitHub issue <https://github.com/castcollab/tesserae2/issues>`_.

Development
-----------

::

    python3 -mvenv venv
    . venv/bin/activate
    pip install -e .

Tests
`````

::

    # run all tests
    tox

    # run only (fast) unit tests
    tox -- -k unit

Building the docs
`````````````````

