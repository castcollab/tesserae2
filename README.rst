Cortexpy
========

.. image:: https://travis-ci.org/winni2k/cortexpy.svg?branch=master
  :target: https://travis-ci.org/winni2k/cortexpy
.. image:: https://coveralls.io/repos/github/winni2k/cortexpy/badge.svg?branch=admin_category
  :target: https://coveralls.io/github/winni2k/cortexpy?branch=admin_category
.. image:: https://landscape.io/github/winni2k/cortexpy/master/landscape.svg?style=flat
  :target: https://landscape.io/github/winni2k/cortexpy/master

Supporting Python versions 3.5 and 3.6.

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

