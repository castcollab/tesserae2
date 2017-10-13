import attr
import pytest

import pycortex.graph as graph
import pycortex.test.builder as builder
import pycortex.graph.serializer as serializer
from pycortex.graph.serializer import collapse_kmer_unitigs


class TestCollapseKmerUnitigsCreatesSingleUnitig(object):
    def test_with_missing_kmer(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3)
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTT')

        # when
        expect = CollapsedKmerUnitgGraphExpectation(collapse_kmer_unitigs(kmer_graph))

        # then
        assert expect.has_kmer_names('GTT')

    def test_with_one_kmer(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('GTT')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTT')

        # when
        unitig_collapsed_kmer_graph = collapse_kmer_unitigs(kmer_graph)

        # then
        assert set(unitig_collapsed_kmer_graph) == {'GTT'}

    def test_with_two_unlinked_kmers(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3).with_kmer('GTT').with_kmer('TTG')
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTG')

        # when
        unitig_collapsed_kmer_graph = collapse_kmer_unitigs(kmer_graph)

        # then
        assert set(unitig_collapsed_kmer_graph) == {'GTT', 'TTG'}

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
        assert expect.has_kmer_names('AAAC')

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
        assert expect.has_kmer_names('AAACC')


@attr.s(slots=True)
class CollapsedKmerUnitgGraphExpectation(object):
    graph = attr.ib()

    def has_kmer_names(self, *kmer_names):
        assert {data['name'] for _, data in self.graph.nodes.data()} == set(kmer_names)
        return self


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
        expect.has_kmer_names('AAAC', 'GAC')

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
        expect.has_kmer_names('AACC', 'C', 'GA')


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
        kmer_graph = serializer.make_graph_json_representable(kmer_graph)

        # then
        assert kmer_graph.node['TTT']['coverage'] == [1, 1]
        assert kmer_graph.node['GTT']['coverage'] == [1, 0]
        assert kmer_graph.node['ATT']['is_missing']

        ttt_node = kmer_graph.node['TTT']
        with pytest.raises(KeyError):
            ttt_node['kmer']
