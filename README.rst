========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis|
        | |coveralls| |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/cortexpy/badge/?style=flat
    :target: https://readthedocs.org/projects/cortexpy
    :alt: Documentation Status

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

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/cortexpy/0.40.0.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/cortexpy/compare/0.40.0...master

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

Set up the project::

    make init

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox

Bugs
====

Please raise a github issue for any bugs.

