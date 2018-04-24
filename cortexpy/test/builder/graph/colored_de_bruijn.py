from unittest.mock import Mock

import attr
from delegation import SingleDelegated

from cortexpy.graph import colored_de_bruijn
from cortexpy.graph.parser.kmer import EmptyKmerBuilder
from cortexpy.graph.parser.streaming import load_de_bruijn_graph
from cortexpy.test.builder import Graph


@attr.s(slots=True)
class ColoredDeBruijnGraphBuilder(object):
    graph = attr.ib(attr.Factory(colored_de_bruijn.ColoredDeBruijnDiGraph))
    colors = attr.ib(init=False)
    kmer_builder = attr.ib(attr.Factory(EmptyKmerBuilder))

    def __attrs_post_init__(self):
        self.with_colors(0)

    def with_node_coverage(self, node, coverage):
        if isinstance(coverage, int):
            coverage = (coverage,)
        else:
            coverage = tuple(coverage)
        assert len(self.colors) == len(coverage)
        self.graph.node[node].coverage = coverage
        return self

    def with_node_kmer(self, node, kmer):
        self.graph.add_node(node, kmer=kmer)
        return self

    def with_node(self, node):
        return self.add_node(node)

    def add_node(self, node):
        if node not in self.graph:
            self.with_node_kmer(node, self.kmer_builder.build_or_get(node))
        return self

    def with_colors(self, *colors):
        assert len(self.graph) == 0
        self.colors = set(colors)
        self.kmer_builder.num_colors = len(self.colors)
        return self

    def with_color(self, color):
        self.colors.add(color)
        self.kmer_builder.num_colors = len(self.colors)

    def add_edge(self, u, v, color=0, key=None):
        if key is not None:
            color = key
        self.add_edge_with_color(u, v, color)
        return self

    def add_edge_with_color(self, u, v, color):
        assert color in self.colors
        self.add_node(u)
        self.add_node(v)
        self.graph.add_edge(u, v, key=color)
        return self

    def add_path(self, *k_strings, color=0, coverage=0):
        if len(k_strings) == 1 and isinstance(k_strings[0], list):
            k_strings = k_strings[0]
        kmer = self.kmer_builder.build_or_get(k_strings[0])
        for cov_color in range(kmer.num_colors):
            kmer.coverage[cov_color] = coverage
        self.graph.add_node(k_strings[0], kmer=kmer)
        if len(k_strings) > 1:
            for k_string1, k_string2 in zip(k_strings[:-1], k_strings[1:]):
                kmer = self.kmer_builder.build_or_get(k_string2)
                for cov_color in range(kmer.num_colors):
                    kmer.coverage[cov_color] = coverage
                self.graph.add_node(k_string2, kmer=kmer)
                self.add_edge_with_color(k_string1, k_string2, color)

    def build(self):
        add_kmers_to_graph(self.graph, num_colors=len(self.colors))
        return self.graph


def add_kmers_to_graph(graph, *, num_colors=1):
    for node, node_data in graph.nodes(data=True):
        if 'kmer' not in node_data:
            kmer_mock = Mock()
            kmer_mock.coverage = tuple(1 for _ in range(num_colors))
            kmer_mock.kmer = node
            graph.node[node]['kmer'] = kmer_mock


class CdbBuilder(SingleDelegated):

    def build(self):
        return load_de_bruijn_graph(self.delegate.build())


def get_cdb_builder():
    return CdbBuilder(Graph())
