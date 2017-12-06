import attr
import networkx as nx


@attr.s(slots=True)
class Builder(object):
    graph = attr.ib(attr.Factory(nx.MultiDiGraph))
    colors_to_link = attr.ib(None)

    def add_kmer(self, kmer, kmer_string):
        self.graph.add_node(kmer_string, kmer=kmer)
        for color, edge_set in enumerate(kmer.edges):
            if self.colors_to_link is not None and color not in self.colors_to_link:
                continue
            for outgoing_kmer in edge_set.get_outgoing_kmer_strings(kmer_string):
                self.graph.add_edge(kmer_string, outgoing_kmer, key=color)
            for incoming_kmer in edge_set.get_incoming_kmer_strings(kmer_string):
                self.graph.add_edge(incoming_kmer, kmer_string, key=color)
        return self
