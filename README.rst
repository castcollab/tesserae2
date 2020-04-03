.. image:: https://github.com/castcollab/tesserae2/workflows/Tests/badge.svg
    :target: Tests
    :alt: Tests

Overview of Tesserae2
=====================

Tesserae2: Fast recombination-aware global and local alignment.

Current version: 1.1.0

Audience
--------


Free software
-------------

Tesserae2 is free software; you can redistribute it and/or modify it under the
terms of the Apache License version 2.0.  Contributions are welcome. Please join us on
`GitHub <https://github.com/castcollab/tesserae2>`_.


Installation
------------

::

    pip install .


Documentation
-------------

CLI
```

After installing Tesserae, the command-line interface can be called by the following command:

::

   tesserae

Bugs
----

For bugs, please raise a `GitHub issue <https://github.com/castcollab/tesserae2/issues>`_.

Development
-----------

::

    python3 -mvenv venv
    . venv/bin/activate
    pip install -r dev-requirements.txt
    pip install -e .

Linting files
`````````````

::

    # reformat all files in src and tests
    black src tests

    # check pep8 against all files in src and tests
    flake8 src tests

    # lint python code for common errors and codestyle issues
    pylint src

Tests
`````

::

    # run all tests
    # right now that includes pep8, black, and pylint
    tox

    # run only (fast) unit tests
    tox -e unit

    # run only linting
    tox -e lint

Bumping versions
````````````````

When making code changes, please bump the version using `bumpversion`. Please make
patch version increments (`bumpversion patch`) for bug fixes, and minor version
increments (`bumpversion minor`) for feature additions **and** incompatible API changes.

Building the docs
`````````````````

Changes
-------

v1.1.0
``````
- tesserae now uses numpy vectorization to speed up the recurrence calculation
  -- @karljohanw
- tesserae now uses multithreading to parallelize the recurrence calculation across
  targets if you are running python v3.8 -- @karljohanw

