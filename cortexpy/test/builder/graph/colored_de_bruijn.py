from unittest.mock import Mock

import attr
from cortexpy.graph import colored_de_bruijn
from cortexpy.graph.parser.kmer import EmptyKmerBuilder


@attr.s(slots=True)
class ColoredDeBruijnGraphBuilder(object):
    graph = attr.ib(attr.Factory(colored_de_bruijn.ColoredDeBruijn))
    colors = attr.ib(attr.Factory(lambda: {0}))
    kmer_builder = attr.ib(attr.Factory(EmptyKmerBuilder))

    def with_node_coverage(self, node, coverage):
        if isinstance(coverage, int):
            coverage = (coverage,)
        else:
            coverage = tuple(coverage)
        node_dict = self.graph.node[node]
        if 'kmer' not in node_dict:
            node_dict['kmer'] = Mock()
        node_dict['kmer'].coverage = coverage
        return self

    def with_node_kmer(self, node, kmer):
        self.graph.add_node(node, kmer=kmer)
        return self

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

    def add_path(self, iterable, color=0, coverage=0):
        k_strings = list(iterable)
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
