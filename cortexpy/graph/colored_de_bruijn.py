from collections import Collection, Mapping

import copy
from itertools import chain

from delegation import SingleDelegated
import attr

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import connect_kmers
from cortexpy.utils import lexlo


@attr.s(slots=True)
class ColoredDeBruijn(Collection):
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

    def __getitem__(self, item):
        return {n for n in chain(self.pred[item], self.succ[item])}

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
            for out_node in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                        is_lexlo=is_lexlo):
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

    @property
    def succ(self):
        return MultiAdjacencyView(self._nodes, EdgeTraversalOrientation.original)

    @property
    def pred(self):
        return MultiAdjacencyView(self._nodes, EdgeTraversalOrientation.reverse)

    def subgraph(self, kmer_strings):
        """Return a subgraph from kmer_strings"""
        dict_view = DictView(self._nodes, set(kmer_strings))
        return ColoredDeBruijn(dict_view, self.graph)

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return ColoredDeBruijn(copy.copy(self._nodes), self.graph.copy())

    def remove_node(self, node):
        for neighbor in self[node]:
            self.in_edges(neighbor, keys=True)
        del self._nodes[node]


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

    def __getitem__(self, item):
        return self._nodes[item]


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
                        ret = [node, kmer_string]
                        if keys:
                            ret.append(color)
                        if data:
                            ret.append({})
                        yield tuple(ret)


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
                edge_kmers.add(lexlo(out_node))
        return iter(edge_kmers)


@attr.s(slots=True)
class MultiAdjacencyView(object):
    _nodes = attr.ib()
    orientation = attr.ib()

    def __getitem__(self, item):
        kmer = self._nodes[item]
        is_lexlo = kmer.kmer == item
        return AdjancencyView(kmer,
                              kmer_string_is_lexlo=is_lexlo,
                              orientation=self.orientation)


@attr.s(slots=True)
class DictView(Mapping):
    base_dict = attr.ib()
    allowed_keys = attr.ib(None)

    def __getitem__(self, item):
        if item in self.allowed_keys:
            return self.base_dict[item]

    def __iter__(self):
        for key, val in self.base_dict.items():
            if key in self.allowed_keys:
                yield key, val

    def __len__(self):
        num_overlap = len(self.allowed_keys & self.base_dict.keys())
        return len(self.base_dict) - num_overlap

    def __copy__(self):
        return {k: v for k, v in self}
