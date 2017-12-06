import attr
import networkx as nx

from cortexpy.graph.parser.kmer import EmptyKmerBuilder


@attr.s(slots=True)
class Builder(object):
    graph = attr.ib(attr.Factory(nx.MultiDiGraph))

    def add_kmer(self, kmer, kmer_string):
        self.graph.add_node(kmer_string, kmer=kmer)
        for color, edge_set in enumerate(kmer.edges):
            for outgoing_kmer in edge_set.get_outgoing_kmer_strings(kmer_string):
                self.graph.add_node(outgoing_kmer)
                self.graph.add_edge(kmer_string, outgoing_kmer, key=color)
            for incoming_kmer in edge_set.get_incoming_kmer_strings(kmer_string):
                self.graph.add_node(incoming_kmer)
                self.graph.add_edge(incoming_kmer, kmer_string, key=color)
        return self
