import json

import pytest

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
        expect.has_n_kmers(1)
        expect.has_one_kmer_with_repr('GTT').is_missing()

    def test_with_one_kmer(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('AAC', 1)
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(1)
        expect.has_one_kmer_with_repr('GTT').is_not_missing()

    @pytest.mark.xfail(reason='Requires MultiDiGraph implementation')
    def test_with_two_unlinked_kmers(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('CAA').with_kmer('AAC')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTG')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(2)
        expect.has_one_kmer_with_repr('GTT').is_not_missing()
        expect.has_one_kmer_with_repr('G').is_not_missing()
        expect.has_n_edges(1).has_n_missing_edges(1)

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
        expect.has_one_kmer_with_repr('AAACC')


class TestCollapseKmerUnitigs(object):
    def test_with_two_node_path_and_three_node_cycle_results_in_two_unitigs(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAA', 1, '.....C..')
                         .with_kmer('AAC', 1, 'a.....G.')
                         .with_kmer('ACG', 1, 'a.g.A...')
                         .with_kmer('CGA', 1, 'a....C..')
                         .with_kmer('GAC', 1, '.c....G.')
                         )

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AAACGAC')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_kmers('AAAC', 'GAC')
        expect.has_one_kmer_with_repr('AAAC').is_not_missing()
        expect.has_one_kmer_with_repr('GAC').is_not_missing()
        expect.has_n_edges(2)
        expect.has_n_missing_edges(0)

    def test_four_node_path_with_one_node_bubble_in_three_nodes(self):
        # given
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAC', 1, '.....C..')
                         .with_kmer('ACC', 1, 'a....CG.')
                         .with_kmer('CCC', 1, 'a.....G.')
                         .with_kmer('CCG', 1, 'ac..A...')
                         .with_kmer('CGA', 1, '.c......'))

        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('AACCCCGA')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(3)
        for kmer in ['AACC', 'C', 'GA']:
            expect.has_one_kmer_with_repr(kmer).is_not_missing()

    def test_two_kmers_one_kmer_apart_do_not_collapse(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('ACC', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GGTTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(3)
        expect.has_n_missing_kmers(1)
        expect.has_n_kmers_with_repr('T', 2)
        expect.has_one_kmer_with_repr('GGT').is_not_missing()

        expect.has_n_edges(2)
        # expect.has_n_missing_edges(2)

    def test_two_linked_kmers_with_incoming_edge_and_missing_kmer_returns_three_unitigs(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.c...C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTTA')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        expect.has_n_kmers(3)

        # fixme: this should be 1 if using MultiDiGraph
        expect.has_n_missing_kmers(2)
        expect.has_one_kmer_with_repr('A').is_missing()
        # fixme: should be not missing if using MultiDiGraph
        expect.has_one_kmer_with_repr('G').is_missing()
        expect.has_one_kmer_with_repr('GTTT').is_not_missing()

        expect.has_n_edges(2)
        # fixme: reintroduce whith MultiDiGraph expect.has_n_missing_edges(2)


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
            serializer.Serializer(kmer_graph).to_json_serializable()
        )

        # then
        expect.has_n_kmers(1).has_one_kmer_with_repr('GTTT').is_not_missing()

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
            serializer.Serializer(kmer_graph).to_json_serializable()
        )

        # then
        expect.has_n_kmers(3)
        expect.has_n_missing_kmers(1)
        expect.has_one_kmer_with_repr('GGT').is_not_missing()
        expect.has_n_kmers_with_repr('T', 2)

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
            serializer.Serializer(kmer_graph).to_json_serializable()
        )

        # then
        expect.has_n_kmers(2)
        expect.has_one_kmer_with_repr('GTTT').is_not_missing()
        expect.has_one_kmer_with_repr('AA').is_missing()
        expect.has_n_edges(1)
        expect.has_n_missing_edges(1)

    def test_linked_kmers_with_incoming_edge_surrounded_by_missing_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a....C..')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('CAAGTTTAGG')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(
            serializer.Serializer(kmer_graph).to_json_serializable()
        )

        # then
        expect.has_n_kmers(4)
        expect.has_one_kmer_with_repr('TT').is_not_missing()
        expect.has_one_kmer_with_repr('CAAGT').is_missing()
        # fixme: should be not missing if using MultiDiGraph
        expect.has_one_kmer_with_repr('GGT').is_missing()
        expect.has_one_kmer_with_repr('AGG').is_missing()

        expect.has_n_edges(3)
        # fixme: should be 2 if using MultiDiGraph
        # expect.has_n_missing_edges(2)
