========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - tests
      - | |travis|
        | |coveralls| |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|

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

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/cortexpy/0.46.0.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/cortexpy/compare/0.46.0...master

.. |wheel| image:: https://img.shields.io/pypi/wheel/cortexpy.svg
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/cortexpy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/cortexpy.svg
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/cortexpy

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/cortexpy.svg
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/cortexpy


.. end-badges

Cortexpy is a Python package for sequence analysis using linked and colored De Bruijn graphs such as
the ones created by `Mccortex <https://github.com/mcveanlab/mccortex>`_.
This project aims to mirror many of the features contained in
`CortexJDK <https://github.com/mcveanlab/CortexJDK>`_.

* Free software: Apache Software License 2.0

Installation
============

::

    pip install cortexpy


Development
===========

1. Install `conda <https://docs.conda.io/en/latest/miniconda.html>`_.
2. Download development and testing tools::

    conda env create -f environment.lock.yml -n my-dev-environment

3. Activate development environment::

    conda activate my-dev-environment

All remaining commands in the development section need to be run in an activated
conda dev environment.

Tests
~~~~~

::

    make test

Update the dev environment
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    # Create a new env from the high-level requirements file
    conda env create -f environment.yml -n another-dev-env

    # activate the new environment
    conda activate another-dev-env

    # save new env to environment.lock.yml
    make lock

Deploy new cortexpy version to pypi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires access credentials for pypi.

::

    make deploy

Bugs
====

Please raise a github issue for any bugs.

