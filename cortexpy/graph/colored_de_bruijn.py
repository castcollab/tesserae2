from collections import Collection, Mapping
from delegation import SingleDelegated
import attr

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import connect_kmers


class SubscriptableKmer(SingleDelegated, Mapping):

    def __getitem__(self, item):
        assert item == 'kmer'
        return self.delegate

    def __str__(self):
        return str(self.delegate)

    def __iter__(self):
        yield 'kmer'

    def __len__(self):
        return 1


@attr.s(slots=True)
class NodeView(Collection):
    _nodes = attr.ib()

    def __call__(self, *args, data=False, **kwargs):
        if data:
            yield from self._nodes.items()
        else:
            yield from self._nodes.keys()

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def __contains__(self, item):
        return item in self._nodes


@attr.s(slots=True)
class EdgeView(object):
    _nodes = attr.ib()

    def __call__(self, *args, data=False, keys=False, **kwargs):
        yield from self._edge_iter(data=data, keys=keys)

    def __len__(self):
        return len(list(self._edge_iter()))

    def __iter__(self):
        return self._edge_iter()

    def __contains__(self, item):
        return item in self._edge_iter()

    def _edge_iter(self, data=False, keys=False):
        for node, kmer in self._nodes.items():
            is_lexlo = node == kmer.kmer
            for color in kmer.colors:
                for kmer_string in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                               is_lexlo=is_lexlo):
                    if kmer_string in self._nodes:
                        if keys:
                            yield (node, kmer_string, color)
                        else:
                            yield (node, kmer_string)


@attr.s(slots=True)
class AdjancencyView(object):
    kmer = attr.ib()
    kmer_string_is_lexlo = attr.ib()
    orientation = attr.ib()

    def __iter__(self):
        edge_kmers = set()
        for color in self.kmer.colors:
            if self.orientation == EdgeTraversalOrientation.original:
                node_iter = self.kmer.edges[color].get_outgoing_kmer_strings(self.kmer.kmer,
                                                                             is_lexlo=self.kmer_string_is_lexlo)
            else:
                node_iter = self.kmer.edges[color].get_incoming_kmer_strings(self.kmer.kmer,
                                                                             is_lexlo=self.kmer_string_is_lexlo)
            for out_node in node_iter:
                edge_kmers.add(out_node)
        return iter(edge_kmers)


@attr.s(slots=True)
class ColoredBeBruijn(Collection):
    """Stores cortex k-mers and conforms to parts of the interface of networkx.MultiDiGraph"""
    _nodes = attr.ib(attr.Factory(dict))
    graph = attr.ib(attr.Factory(dict))

    @property
    def nodes(self):
        return NodeView(self._nodes)

    @property
    def node(self):
        return self._nodes

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def __contains__(self, item):
        return item in self._nodes

    def add_node(self, kmer_string, *, kmer):
        self._nodes[kmer_string] = SubscriptableKmer(kmer)

    def add_nodes_from(self, node_iterable):
        for node in node_iterable:
            self.add_node(node[0], kmer=node[1])

    def is_multigraph(self):
        return True

    def is_directed(self):
        return True

    @property
    def edges(self):
        return EdgeView(self._nodes)

    def out_edges(self, node, keys=False, default=None):
        for out_kmer in self.succ(node)
            for color in kmer.colors:
            for out_node in kmer.edges[color].get_outgoing_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (node, out_node, color)
                else:
                    yield (node, out_node)

    def in_edges(self, node, keys=False, default=None):
        kmer = self._nodes[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for in_node in kmer.edges[color].get_incoming_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (in_node, node, color)
                else:
                    yield (in_node, node)

    def out_degree(self, node):
        return len(list(self.out_edges(node)))

    def in_degree(self, node):
        return len(list(self.in_edges(node)))

    def add_edge(self, first, second, *, key):
        """Note: edges can only be added to existing nodes"""
        first_kmer = self._nodes[first]
        second_kmer = self._nodes[second]
        connect_kmers(first_kmer, second_kmer, color=key)

    def add_edges_from(self, edge_iterable):
        for edge in edge_iterable:
            self.add_edge(edge[0], edge[1], key=edge[2])

    def __str__(self):
        return '\n'.join(self._nodes.items())

    def succ(self, node):
        kmer = self._nodes[node]
        is_lexlo = kmer.kmer == node
        return AdjancencyView(kmer, kmer_string_is_lexlo=is_lexlo,
                              orientation=EdgeTraversalOrientation.original)

    def pred(self, node):
        kmer = self._nodes[node]

        is_lexlo = kmer.kmer == node
        return AdjancencyView(kmer, kmer_string_is_lexlo=is_lexlo,
                              orientation=EdgeTraversalOrientation.reverse)
