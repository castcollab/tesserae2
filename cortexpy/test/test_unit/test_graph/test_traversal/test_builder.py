import cortexpy.graph as graph
from cortexpy.graph.parser.kmer import EmptyKmerBuilder, connect_kmers
from cortexpy.test.expectation import KmerGraphExpectation


class Test(object):
    def test_without_kmers_builds_empty_graph(self):
        # given
        builder = graph.Builder()

        # when/then
        assert len(builder.graph) == 0

    def test_with_one_kmer_adds_kmer_to_graph(self):
        # given
        builder = graph.Builder()
        ek_builder = EmptyKmerBuilder()

        # when
        expect = KmerGraphExpectation(builder.add_kmer(ek_builder.build_or_get('AAA'), 'AAA').graph)

        # then
        expect.has_n_nodes(1)
        expect.has_node('AAA')
        expect.has_n_edges(0)

    def test_with_two_kmers_adds_edge_to_graph_for_one_color(self):
        # given
        builder = graph.Builder()
        ek_builder = EmptyKmerBuilder(num_colors=1)

        kmer1 = ek_builder.build_or_get('AAA')
        kmer2 = ek_builder.build_or_get('AAT')
        connect_kmers(kmer1, kmer2, color=0)

        builder.add_kmer(kmer1, 'AAA')
        builder.add_kmer(kmer2, 'AAT')

        # when
        expect = KmerGraphExpectation(builder.graph)

        # then
        expect.has_n_nodes(2)
        expect.has_node('AAA').has_coverages(0)
        expect.has_node('AAT').has_coverages(0)
        expect.has_n_edges(1)

    def test_with_two_kmers_adds_edge_to_graph_for_two_colors(self):
        # given
        builder = graph.Builder()
        ek_builder = EmptyKmerBuilder(num_colors=2)

        kmer1 = ek_builder.build_or_get('AAA')
        kmer2 = ek_builder.build_or_get('AAT')
        connect_kmers(kmer1, kmer2, color=1)

        builder.add_kmer(kmer1, 'AAA')
        builder.add_kmer(kmer2, 'AAT')

        # when
        expect = KmerGraphExpectation(builder.graph)

        # then
        expect.has_node('AAA').has_coverages(0)
        expect.has_node('AAT').has_coverages(0)
        expect.has_n_nodes(2)
        expect.has_edge('AAA', 'AAT', 1)
        expect.has_n_edges(1)

    def test_with_two_revcomp_kmers_adds_edge_to_graph_for_one_color(self):
        # given
        builder = graph.Builder()
        ek_builder = EmptyKmerBuilder(num_colors=1)

        kmer1 = ek_builder.build_or_get('AAA')
        kmer2 = ek_builder.build_or_get('AAT')
        connect_kmers(kmer1, kmer2, color=0)

        builder.add_kmer(kmer1, 'TTT')
        builder.add_kmer(kmer2, 'ATT')

        # when
        expect = KmerGraphExpectation(builder.graph)

        # then
        expect.has_nodes('TTT', 'ATT')
        expect.has_n_edges(1)
