from unittest.mock import Mock

import attr
import networkx as nx


@attr.s()
class NetworkxGraphBuilder(object):
    graph = attr.ib(attr.Factory(nx.DiGraph))

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

    def build(self):
        self.graph = add_kmers_to_graph(self.graph)
        return self.graph


def add_kmers_to_graph(graph):
    out_graph = graph.copy()
    for node, node_data in graph.nodes(data=True):
        if 'kmer' not in node_data:
            kmer_mock = Mock()
            kmer_mock.coverage = (1,)
            kmer_mock.kmer = node
            out_graph.node[node]['kmer'] = kmer_mock
    return out_graph
