from collections import Collection, Mapping

import copy
from itertools import chain

from delegation import SingleDelegated
import attr

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import connect_kmers
from cortexpy.utils import lexlo


def build_cdb_graph_from_header(header, kmer_generator=None):
    return build_cdb_graph(sample_names=header.sample_names,
                           kmer_size=header.kmer_size,
                           num_colors=header.num_colors,
                           colors=header.colors,
                           kmer_generator=kmer_generator)


def build_cdb_graph_from_ra_parser(ra_parser, kmer_generator=None):
    return build_cdb_graph(sample_names=ra_parser.sample_names,
                           kmer_size=ra_parser.kmer_size,
                           num_colors=ra_parser.num_colors,
                           colors=ra_parser.colors)


def build_cdb_graph(*, sample_names, kmer_size, num_colors, colors, kmer_generator=None):
    """Colored de Bruijn graph constructor"""
    if kmer_generator is not None:
        graph = ColoredDeBruijnDiGraph({k.kmer: SubscriptableKmer(k) for k in kmer_generator})
    else:
        graph = ColoredDeBruijnDiGraph()
    graph.graph['sample_names'] = sample_names
    graph.graph['kmer_size'] = kmer_size
    graph.graph['num_colors'] = num_colors
    graph.graph['colors'] = colors
    return graph


@attr.s(slots=True)
class ColoredDeBruijnDiGraph(Collection):
    """Stores cortex k-mers and conforms to parts of the interface of networkx.MultiDiGraph"""
    _nodes = attr.ib(attr.Factory(dict))
    graph = attr.ib(attr.Factory(dict))

    @property
    def nodes(self):
        return NodeView(self._nodes)

    @property
    def node(self):
        return NodeView(self._nodes)

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def __contains__(self, item):
        return item in self._nodes.keys()

    def __getitem__(self, item):
        return self.succ[item]

    def add_node(self, kmer_string, *, kmer):
        self._nodes[kmer_string] = SubscriptableKmer(kmer)

    def add_nodes_from(self, node_iterable):
        for node in node_iterable:
            self.add_node(node[0], kmer=node[1])

    def is_multigraph(self):
        return True

    def is_directed(self):
        return True

    def is_conistent(self):
        return False

    @property
    def edges(self):
        return EdgeView(self)

    def out_edges(self, node, keys=False, default=None, data=None):
        kmer = self.node[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for out_node in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                        is_lexlo=is_lexlo):
                if keys:
                    yield (node, out_node, color)
                else:
                    yield (node, out_node)

    def in_edges(self, node, keys=False, default=None, data=None):
        kmer = self.node[node]
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
        first_kmer = self.node[first]
        second_kmer = self.node[second]
        connect_kmers(first_kmer, second_kmer, color=key)

    def add_edges_from(self, edge_iterable):
        for edge in edge_iterable:
            self.add_edge(edge[0], edge[1], key=edge[2])

    def __str__(self):
        return '\n'.join(self._nodes.keys())

    @property
    def succ(self):
        return MultiAdjacencyView(self._nodes, EdgeTraversalOrientation.original)

    @property
    def pred(self):
        return MultiAdjacencyView(self._nodes, EdgeTraversalOrientation.reverse)

    def subgraph(self, kmer_strings):
        """Return a subgraph from kmer_strings"""
        dict_view = DictView(self._nodes, set(kmer_strings))
        return ColoredDeBruijnDiGraph(dict_view, self.graph)

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return ColoredDeBruijnDiGraph(copy.copy(self._nodes), self.graph.copy())

    def remove_node(self, node):
        try:
            node_kmer = self.node[node]
        except KeyError:
            return
        for succ in self.succ[node]:
            succ_kmer = self.node[succ]
            letter, _ = succ_kmer.find_letters_of_edge_from_kmer(node_kmer)
            for color in succ_kmer.colors:
                succ_kmer.edges[color].remove_edge(letter)
        for pred in self.pred[node]:
            pred_kmer = self.node[pred]
            letter, _ = pred_kmer.find_letters_of_edge_to_kmer(node_kmer)
            for color in pred_kmer.colors:
                pred_kmer.edges[color].remove_edge(letter)
        del self._nodes[node]

    def remove_nodes_from(self, nodes):
        for node in nodes:
            self.remove_node(node)


class ColoredDeBruijnGraph(SingleDelegated):

    def is_directed(self):
        return False

    def __getitem__(self, item):
        return {
            n: c for n, c in
            chain(((k, self.delegate.pred[item][k]) for k in self.delegate.pred[item]),
                  ((k, self.delegate.succ[item][k]) for k in self.delegate.succ[item]))
        }

    def __len__(self):
        """Wish this was picked up by SingleDelegated ..."""
        return len(self.delegate)


@attr.s(slots=True)
class ColoredDeBruijnLexloDiGraph(SingleDelegated):

    def add_node(self, kmer_string, *, kmer):
        kmer_string = lexlo(kmer_string)
        self.delegate.node[kmer_string] = SubscriptableKmer(kmer)


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


class LexloNodeView(SingleDelegated):

    def __getitem__(self, item):
        return self._nodes[lexlo(item)]


@attr.s(slots=True)
class EdgeView(object):
    graph = attr.ib()

    def __call__(self, *args, data=False, keys=False, **kwargs):
        node = None
        if args:
            node = args[0]
        if not node:
            yield from self._edge_iter(data=data, keys=keys)
        else:
            if self.graph.is_directed():
                yield from self.graph.out_edges(node, data=data, keys=keys)
            else:
                yield from self._edge_iter(data=data, keys=keys)

    def __len__(self):
        return len(list(self._edge_iter()))

    def __iter__(self):
        return self._edge_iter()

    def __contains__(self, item):
        return item in self._edge_iter()

    def _edge_iter(self, data=False, keys=False):
        for node, kmer in self.graph.nodes(data=True):
            is_lexlo = node == kmer.kmer
            for color in kmer.colors:
                for kmer_string in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                               is_lexlo=is_lexlo):
                    if kmer_string in self.graph.nodes:
                        ret = [node, kmer_string]
                        if keys:
                            ret.append(color)
                        if data:
                            ret.append({})
                        yield tuple(ret)


@attr.s(slots=True)
class AdjancencyView(object):
    _nodes = attr.ib()
    kmer = attr.ib()
    kmer_string_is_lexlo = attr.ib()
    orientation = attr.ib()

    @property
    def colors(self):
        return self.kmer.colors

    def __iter__(self):
        edge_kmers = set()
        for color in self.colors:
            if self.orientation == EdgeTraversalOrientation.original:
                node_iter = self.kmer.edges[color].get_outgoing_kmer_strings(
                    self.kmer.kmer, is_lexlo=self.kmer_string_is_lexlo
                )
            else:
                node_iter = self.kmer.edges[color].get_incoming_kmer_strings(
                    self.kmer.kmer, is_lexlo=self.kmer_string_is_lexlo
                )
            for out_node in node_iter:
                edge_kmers.add(lexlo(out_node))
        return iter(edge_kmers)

    def __getitem__(self, item):
        other = self._nodes[item]
        edge_colors = set()
        if self.orientation == EdgeTraversalOrientation.original:
            edge_func = self.kmer.has_outgoing_edge_to_kmer_in_color
        else:
            edge_func = self.kmer.has_incoming_edge_from_kmer_in_color
        try:
            for color in self.colors:
                if edge_func(other, color):
                    edge_colors.add(color)
        except ValueError:
            pass
        return edge_colors


@attr.s(slots=True)
class MultiAdjacencyView(object):
    _nodes = attr.ib()
    orientation = attr.ib()

    def __getitem__(self, item):
        kmer = self._nodes[item]
        is_lexlo = kmer.kmer == item
        return AdjancencyView(self._nodes,
                              kmer,
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
