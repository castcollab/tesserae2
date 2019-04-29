On link-informed graph traversal
==================================

The all-simple-paths algorithm as implemented in :py:func:`nx.all_simple_paths` inside the
:py:mod:`networkx` (abbreviated below as :py:mod:`nx`) package version 2.2 uses a
depth-first traversal scheme to find all possible paths from a start node to one or more
end nodes [#]_.

For example, let nodes ``A``-``F`` represent unitigs in a De Bruijn graph created from sequencing
reads of transcripts:

.. mermaid::

    graph LR;
        A-->B
        A-->C
        B-->D
        C-->D
        D-->E
        D-->F

An all-simple-paths traversal starting at ``A`` will return the paths ``ABDE``, ``ABDF``, ``ACDE``,
and ``ACDF``. However, what if the sequenced reads that were used to create this graph only
originated from two paths: ``ABDE`` and
``ACDF``? Can some of these sequencing reads be used to restrict the paths returned by the
all-simple-paths algorithm?

Mccortex_ provides a data structure called "links_" for annotating De Bruijn graphs in Cortex_
format. In the example above, links can be used to store information on a read that covers both the
``A`` --> ``B`` and ``D`` --> ``E`` junctions. Cortexpy can use these links to performed a
"link-informed" (that is, a restricted) all-simple-paths traversal.

.. _Cortex: https://github.com/iqbal-lab/cortex
.. _Mccortex: https://github.com/mcveanlab/mccortex
.. _links: https://github.com/mcveanlab/mccortex/wiki/Graph-links


.. [#] Cortexpy currently uses a copy of :py:func:`nx.all_simple_paths`
       named :py:func:`~cortexpy.graph.interactor._all_simple_paths_graph`.


Link-informed graph traversal in :py:mod:`cortexpy`
---------------------------------------------------

Cortexpy uses :py:mod:`networkx` algorithms
````````````````````````````````````````````

Cortexpy represents Cortex graphs as :py:class:`nx.DiGraph` objects [#]_. This
allows the easy application of :py:mod:`networkx` algorithms to Cortex graphs.
:py:mod:`cortexpy` achieves link-informed traversal by wrapping a Cortex graph in a
:py:class:`~cortexpy.links.LinkedGraphTraverser` object, which modifies the behavior of the
:py:meth:`~cortexpy.links.LinkedGraphTraverser.__getitem__` method. To understand why this works, let us first take a look at the
all-simple-paths algorithm:

.. code-block:: python3
   :lineno-start: 1
   :emphasize-lines: 22

    def _all_simple_paths_graph(G, source, target, cutoff=None):
        """This function was copied from Networkx before being edited by Warren Kretzschmar
        todo: switch back to nx.all_simple_paths once Networkx 2.2 is released"""
        assert cutoff is not None
        if target in G:
            targets = {target}
        else:
            targets = set(target)
        visited = collections.OrderedDict.fromkeys([source])
        stack = [iter(G[source])]
        while stack:
            children = stack[-1]
            child = next(children, None)
            if child is None:
                stack.pop()
                visited.popitem()
            elif len(visited) < cutoff:
                if child in targets:
                    yield list(visited) + [child]
                elif child not in visited:
                    visited[child] = None
                    stack.append(iter(G[child]))
            else:  # len(visited) == cutoff:
                if child in targets or len(targets & children) != 0:
                    yield list(visited) + [child]
                stack.pop()
                visited.popitem()

The key line here is the highlighted line 22. This is the line that appends an iterator of a node's
successors to the stack of nodes to visit. The algorithm asks the graph object :py:obj:`G`
for the successor nodes of :py:obj:`child` by calling the
:py:meth:`~cortexpy.links.LinkedGraphTraverser.__getitem__` method of :py:obj:`G`::

    G[child]

This means that we can restrict the paths returned by :py:func:`_all_simple_paths_graph` by
restricting the successor nodes returned by :py:obj:`G`.

.. [#] The implementation is not perfect and could use some improvement.

:py:class:`~cortexpy.links.LinkedGraphTraverser` restricts simple paths using links
```````````````````````````````````````````````````````````````````````````````````

The :py:meth:`~cortexpy.links.LinkedGraphTraverser.__getitem__` method of
:py:class:`~cortexpy.links.LinkedGraphTraverser` restricts the returned
successors using the following rules:

1. If no link information exists for the query node, then return all successors.
2. Otherwise, return only successors that are consistent with the links encountered on the path from
   start to query node.
3. If the query node is annotated with links, pick up all links.
4. For each successor node, only retain links that are consistent with the path taken from the start
   to this successor node.
5. For each successor node, drop links that are no longer relevant to the successor node
   (i.e. links that have expired)
