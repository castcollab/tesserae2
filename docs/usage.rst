Installation
============

Pip
---

::

    pip install cortexpy

Conda
-----

::

    conda install cortexpy


Getting started
===============

Inspecting cortex files
-----------------------

After you have installed cortexpy, you can start using it
to access cortex files from python::

    from cortexpy.graph.parser.random_access import RandomAccess

    # make sure to open the cortex graph in binary mode
    with open('my_graph.ctx', 'rb') as fh:
        ra = RandomAccess(fh)

        # let's see if our favorite kmer is in the graph
        if 'AAA' in ra:
            print('AAA exists in my_graph.ctx!')

