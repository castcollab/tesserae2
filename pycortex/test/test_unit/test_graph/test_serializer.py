import json

import pycortex.graph as graph
import pycortex.graph.serializer as serializer
import pycortex.test.builder as builder
from pycortex.graph.serializer import collapse_kmer_unitigs
from pycortex.test.expectation.kmer import CollapsedKmerUnitgGraphExpectation


class TestCollapseKmerUnitigsCreatesSingleUnitig(object):
    def test_with_missing_kmer(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3)
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        assert expect.has_kmers('GTT')

    def test_with_one_kmer(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('GTT')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_kmers('GTT')

    def test_with_two_unlinked_kmers(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('GTT').with_kmer('TTG')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTG')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_kmers('GTT', 'TTG')

    def test_with_two_linked_kmers(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAA', 0, '.....C..')
                         .with_kmer('AAC', 0, 'a.......'))

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        assert expect.has_kmers('AAAC')

    def test_with_three_linked_kmers(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAA', 0, '.....C..')
                         .with_kmer('AAC', 0, 'a....C..')
                         .with_kmer('ACC', 0, 'a.......'))

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AAACC')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(1)
        expect.has_kmer('AAACC')


class TestCollapseKmerUnitigs(object):
    def test_with_two_node_path_and_three_node_cycle_results_in_two_unitigs(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAA', 0, '.....C..')
                         .with_kmer('AAC', 0, 'a.....G.')
                         .with_kmer('ACG', 0, 'a.g.A...')
                         .with_kmer('CGA', 0, 'a....C..')
                         .with_kmer('GAC', 0, '.c....G.')
                         )

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AAACGAC')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_kmers('AAAC', 'GAC')

    def test_four_node_path_with_one_node_bubble_in_three_nodes(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAC', 0, '.....C..')
                         .with_kmer('ACC', 0, 'a....CG.')
                         .with_kmer('CCC', 0, 'a.....G.')
                         .with_kmer('CCG', 0, 'ac..A...')
                         .with_kmer('CGA', 0, '.c......'))

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AACCCCGA')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(3)
        for kmer in ['AACC', 'C', 'GA']:
            expect.has_kmer(kmer)


class TestMakeGraphJsonRepresentable(object):
    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers_jsonifiable(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(
            serializer.make_graph_json_representable(kmer_graph))

        # then
        expect.has_n_kmers(3)
        expect.has_kmer('TTT').has_coverages(1, 1).is_not_missing()
        expect.has_kmer('GTT').has_coverages(1, 0).is_not_missing()
        expect.has_kmer('ATT').is_missing()

        expect.has_n_edges(2)
        expect.has_n_missing_edges(1)


class TestToJson(object):
    def test_two_linked_kmers_are_jsonifiable(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        kmer_json = serializer.Serializer(kmer_graph).to_json()

        # then
        json.loads(kmer_json)  # does not raise


class TestToJsonSerializable(object):
    def test_two_linked_kmers_collapse_to_one_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(
            serializer.Serializer(kmer_graph, collapse_kmer_unitigs=True).to_json_serializable()
        )

        # then
        expect.has_n_kmers(1).has_kmer('GTTT').is_not_missing()

    def test_two_kmers_one_kmer_apart_do_not_collapse(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('ACC', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GGTTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(
            serializer.Serializer(kmer_graph, collapse_kmer_unitigs=True).to_json_serializable()
        )

        # then
        expect.has_n_kmers(3)
        expect.has_kmer('TTT').is_not_missing()
        expect.has_kmer('GTT').is_missing()
        expect.has_kmer('GGT').is_not_missing()

        expect.has_n_edges(2)
        expect.has_n_missing_edges(2)

    def test_unlinked_kmers_followed_by_two_linked_kmers_collapse_to_two_unitigs(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTTAA')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(
            serializer.Serializer(kmer_graph, collapse_kmer_unitigs=True).to_json_serializable()
        )

        # then
        expect.has_n_kmers(2)
        expect.has_kmer('GTTT').is_not_missing()
        expect.has_kmer('TTAA').is_missing()

