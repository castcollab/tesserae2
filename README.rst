Cortexpy
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        | |coveralls| |codecov|
        | |codeclimate|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/cortexpy/badge/?style=flat
    :target: https://readthedocs.org/projects/cortexpy
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/winni2k/cortexpy.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/winni2k/cortexpy

.. |requires| image:: https://requires.io/github/winni2k/cortexpy/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/winni2k/cortexpy/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/winni2k/cortexpy/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/winni2k/cortexpy

.. |codecov| image:: https://codecov.io/github/winni2k/cortexpy/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/winni2k/cortexpy

.. |codeclimate| image:: https://codeclimate.com/github/winni2k/cortexpy/badges/gpa.svg
   :target: https://codeclimate.com/github/winni2k/cortexpy
   :alt: CodeClimate Quality Status

.. |version| image:: https://img.shields.io/pypi/v/cortexpy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/cortexpy

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/cortexpy/v0.32.2.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/cortexpy/compare/v0.32.2...master

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


.. image:: https://travis-ci.org/winni2k/cortexpy.svg?branch=master
  :target: https://travis-ci.org/winni2k/cortexpy
.. image:: https://coveralls.io/repos/github/winni2k/cortexpy/badge.svg?branch=admin_category
  :target: https://coveralls.io/github/winni2k/cortexpy?branch=admin_category
.. image:: https://landscape.io/github/winni2k/cortexpy/master/landscape.svg?style=flat
  :target: https://landscape.io/github/winni2k/cortexpy/master

Cortexpy is a Python package for sequence analysis using linked and colored de Bruijn graphs such as
the ones created by `Mccortex <https://github.com/mcveanlab/mccortex>`_.
This project aims to mirror many of the features contained in
`CortexJDK <https://github.com/mcveanlab/CortexJDK>`_.


Install
-------

Install the latest version of cortexpy::

    pip install cortexpy

Bugs
----

Please raise a github issue for any bugs.

Development
-----------

Set up the project::

    make init

See `cortexpy/__main__.py` for example uses of cortexpy.

Publishing to Pypi
------------------

::

    # assuming version 0.0.1

    # Step 1: Set __version__ in cortexpy/__init__.py

    # Step 2: Make version commit
    git commit -am'Bump to version 0.0.1'

    # Step 3: tag commit
    git tag 0.0.1

    # Step 4: publish to pypi
    make deploy

