Overview of Tesserae2
=====================

Tesserae2: Fast recombination-aware global and local alignment.

Audience
--------


Free software
-------------

Tesserae2 is free software; you can redistribute it and/or modify it under the
terms of the Apache License version 2.0.  Contributions are welcome. Please join us on `GitHub <https://github.com/winni2k/cortexpy>`_.


Installation
------------

::

    pip install tesserae


Documentation
-------------

For more information, please see the tesserae documentation_.

.. _documentation: https://tesserae.readthedocs.io/en/latest/index.html

Bugs
----

For bugs, please raise a `GitHub issue <https://github.com/winni2k/cortexpy/issues>`_.

Development
-----------

::

    python3 -mvenv venv
    . venv/bin/activate
    pip install -e .

Tests
`````

::

    make test

Building the docs
`````````````````

The documentation is automatically built by read-the-docs on push to master.
To build the documentation manually::

    # install sphinx dependencies
    pip install -r docs/requirements.txt

    make docs
