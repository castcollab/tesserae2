Overview of Cortexpy_
=====================

.. start-badges

.. list-table::
    :stub-columns: 1

    * - tests
      - | |travis|
        | |coveralls| |codecov|
    * - package
      - | |version| |wheel| |supported-versions|
        | |supported-implementations| |commits-since|
    * - docs
      - | |readthedocs|

.. |travis| image:: https://travis-ci.org/winni2k/cortexpy.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/winni2k/cortexpy

.. |coveralls| image:: https://coveralls.io/repos/winni2k/cortexpy/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/winni2k/cortexpy

.. |codecov| image:: https://codecov.io/github/winni2k/cortexpy/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/winni2k/cortexpy

.. |version| image:: https://img.shields.io/pypi/v/cortexpy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/cortexpy

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/cortexpy/v0.46.1.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/cortexpy/compare/v0.46.1...master

.. |wheel| image:: https://img.shields.io/pypi/wheel/cortexpy.svg
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/cortexpy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/cortexpy.svg
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/cortexpy

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/cortexpy.svg
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/cortexpy

.. |readthedocs| image:: https://readthedocs.org/projects/cortexpy/badge/?version=latest
   :target: https://cortexpy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. end-badges

Cortexpy is a Python package for sequence analysis using linked and colored De Bruijn graphs such as
the ones created by Cortex_ and Mccortex_.
This project aims to mirror many of the features contained in
`CortexJDK <https://github.com/mcveanlab/CortexJDK>`_.

.. _cortexpy: https://github.com/winni2k/cortexpy

Cortexpy also comes with a command-line tool for basic inspection and manipulation of Cortex graphs with and without links.

.. _Cortex: https://github.com/iqbal-lab/cortex
.. _Mccortex: https://github.com/mcveanlab/mccortex

Audience
--------

The audience of cortexpy is researchers working with colored De Bruijn graphs and link information in Cortex_ and Mccortex_ format.


Free software
-------------

Cortexpy is free software; you can redistribute it and/or modify it under the
terms of the :doc:`Apache License version 2.0</license>`.  Contributions are welcome. Please join us on `GitHub <https://github.com/winni2k/cortexpy>`_.


Installation
------------

::

    pip install cortexpy


Documentation
-------------

For more information, please see cortexpy documentation_.

.. _documentation: https://cortexpy.readthedocs.io/en/latest/index.html

Bugs
----

Please raise a github issue for any bugs.

Development
-----------

1. Install `conda <https://docs.conda.io/en/latest/miniconda.html>`_.
2. Download development and testing tools::

    conda env create -f environment.yml -n my-dev-environment

3. Activate development environment::

    conda activate my-dev-environment

All remaining commands in the development section need to be run in an activated
conda dev environment.



Tests
~~~~~

::

    make test

Deploy new cortexpy version to pypi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires access credentials for pypi.

::

    make deploy

Building the docs
~~~~~~~~~~~~~~~~~

The documentation is automatically built by read-the-docs on push to master.
To build the documentation manually::

    # install sphinx dependencies
    pip install docs/requirements.txt

    make docs

Update the dev environment
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section is experimental because it does not work on travis-CI yet.

::

    # Create a new env from the high-level requirements file
    conda env create -f environment.yml -n another-dev-env

    # activate the new environment
    conda activate another-dev-env

    # save new env to environment.lock.yml
    make lock

