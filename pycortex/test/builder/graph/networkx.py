from unittest.mock import Mock

import attr
import networkx as nx


@attr.s(slots=True)
class NetworkxGraphBuilder(object):
    graph = attr.ib(attr.Factory(nx.MultiDiGraph))
    colors = attr.ib(attr.Factory(set))

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

    def with_colors(self, *colors):
        self.colors = set(colors)
        return self

    def add_edge_with_color(self, u, v, color):
        self.colors.add(color)
        self.graph.add_edge(u, v, key=color)
        return self

    def build(self):
        self.graph = add_kmers_to_graph(self.graph, num_colors=len(self.colors))
        return self.graph


def add_kmers_to_graph(graph, *, num_colors=1):
    out_graph = graph.copy()
    for node, node_data in graph.nodes(data=True):
        if 'kmer' not in node_data:
            kmer_mock = Mock()
            kmer_mock.coverage = tuple(1 for _ in range(num_colors))
            kmer_mock.kmer = node
            out_graph.node[node]['kmer'] = kmer_mock
    return out_graph
