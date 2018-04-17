from collections import Collection, Mapping
from delegation import SingleDelegated
import attr

from cortexpy.graph.parser.kmer import connect_kmers


class SubscriptableKmer(SingleDelegated):

    def __contains__(self, item):
        if item == 'kmer':
            return True
        return False

    def __getitem__(self, item):
        assert item == 'kmer'
        return self

    def __str__(self):
        return str(self.delegate)


@attr.s(slots=True)
class NodeView(Collection):
    _nodes = attr.ib()

    def __call__(self, *args, **kwargs):
        yield from self._nodes.items()

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def __contains__(self, item):
        return item in self._nodes


@attr.s(slots=True)
class EdgeView(object):
    _nodes = attr.ib()

    def __call__(self, *args, **kwargs):
        yield from self._edge_iter()

    def __len__(self):
        return len(list(self._edge_iter()))

    def __iter__(self):
        return self._edge_iter()

    def __contains__(self, item):
        return item in self._edge_iter()

    def _edge_iter(self):
        for node, kmer in self._nodes.items():
            is_lexlo = node == kmer.kmer
            for color in kmer.colors:
                for kmer_string in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                               is_lexlo=is_lexlo):
                    if kmer_string in self._nodes:
                        yield (node, kmer_string, color)


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
        kmer = self._nodes[node]
        is_lexlo = kmer.kmer == node
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
