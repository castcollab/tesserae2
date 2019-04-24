********
Tutorial
********

The cortexpy package consists of a python API and a command-line tool for working with Cortex graphs.
Below, we start by looking at how to use the python API to perform an example workflow.

Using the python API to filter Cortex graphs
============================================

.. _`Building Cortex files`:

Building Cortex files
---------------------

.. _Mccortex: https://github.com/mcveanlab/mccortex
.. _biopython: https://pypi.org/project/biopython/

Let's start by by creating two Cortex files to work with. At present, cortexpy does not provide a way
to easily create a Cortex file, so we will instead use Mccortex_. Mccortex can be compiled from source
or installed using `bioconda <https://bioconda.github.io/recipes/mccortex/README.html>`_.

Let's start by creating two FASTA files from which to create two Cortex files::

    echo -e '>1\nAAAAA' > file1.fasta
    echo -e '>1\nCCCCC' > file2.fasta

We now have two FASTA files each containing a single sequence.  We can now build a Cortex graph
from each file. We choose to use a kmer-size of 5::

    mccortex 5 build --sort -k 5 --sample file1 -1 file1.fasta file1.ctx
    mccortex 5 build --sort -k 5 --sample file2 -1 file2.fasta file2.ctx

We now have two cortex files: ``file1.ctx`` and ``file2.ctx``. As the Cortex format represents
colored De Bruijn graphs, we could have stored the information from the two FASTA files in a single
graph as two separate colors. However, we are creating two files in order to demonstrate the
cortexpy API later on.

We can check what kmers are stored in each graph using the :command:`cortexpy` command-line tool::

    > cortexpy view graph file1.ctx
    AAAAA 1 ........

    > cortexpy view graph file2.ctx
    CCCCC 1 ........

This output tells us that each graph consists of a single kmer with coverage 1.


Inspecting Cortex graphs in Python
----------------------------------

Cortexpy offers many ways to inspect Cortex files. Much of that functionality is available through
the :py:class:`RandomAccess` class. Let us start by loading a Cortex file inside python::

    >>> from cortexpy.graph.parser.random_access import RandomAccess
    >>> # make sure to open the cortex graph in binary mode
    >>> ra = RandomAccess(open('file1.ctx', 'rb'))

We can now interrogate the ``ra`` object. Let's see what the header size of the Cortex file is::

    >>> ra.header.kmer_size
    5

Let's check if the kmer ``AAAAA`` exists in the graph and retrieve it::

    >>> 'AAAAA' in ra
    True
    >>> ra['AAAAA']
    Kmer(_kmer_data=KmerData(_data=b'\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00', kmer_size=5, num_colors=1, _kmer='AAAAA', _coverage=None, _edges=None), num_colors=1, kmer_size=5, _revcomp=None)

We can see that the returned kmer object contains information on the kmer size (5) and the number of colors stored in the kmer (1).

Now let's put it all together and search both graphs that we created while `Building Cortex files`_ for our kmer of interest, ``AAAAA``:

.. code-block:: python
   :caption: search.py
   :name: search-py

    from cortexpy.graph.parser.random_access import RandomAccess

    for graph in ['file1.ctx', 'file2.ctx']:
        # make sure to open the cortex graph in binary mode
        with open(graph, 'rb') as fh:
            ra = RandomAccess(fh)

            # let's see if our favorite kmer is in the graph
            if 'AAAAA' in ra:
                print(f'AAAAA exists in {graph}!')

This is what we see if we run this code from the command line::

    > python3 search.py
    AAAAA exists in file1.ctx!
