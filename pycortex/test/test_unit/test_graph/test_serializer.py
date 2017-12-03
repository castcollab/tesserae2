import json

import attr

import pycortex.graph as graph
import pycortex.graph.serializer as serializer
import pycortex.test.builder as builder
from pycortex.test.expectation.kmer import CollapsedKmerUnitgGraphExpectation


@attr.s(slots=True)
class SerializerTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    contig_to_retrieve = attr.ib(None)
    retriever = attr.ib(None)

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_contig_to_retrieve(self, contig):
        self.contig_to_retrieve = contig
        return self

    def run(self):
        self.retriever = graph.ContigRetriever(self.graph_builder.build())
        return self.retriever.get_kmer_graph(self.contig_to_retrieve)


@attr.s(slots=True)
class CollapseKmerUnitigsTestDriver(object):
    serializer_driver = attr.ib(attr.Factory(SerializerTestDriver))

    def __getattr__(self, name):
        serializer_method = getattr(self.serializer_driver, name)

        def method(*args):
            serializer_method(*args)
            return self

        return method

    def run(self):
        kmer_graph = self.serializer_driver.run()
        collapser = (serializer
                     .UnitigCollapser(kmer_graph, colors=self.serializer_driver.retriever.colors)
                     .collapse_kmer_unitigs())
        return CollapsedKmerUnitgGraphExpectation(collapser.unitig_graph)


@attr.s(slots=True)
class JsonSerializableTestDriver(object):
    serializer_driver = attr.ib(attr.Factory(SerializerTestDriver))

    def __getattr__(self, name):
        serializer_method = getattr(self.serializer_driver, name)

        def method(*args):
            serializer_method(*args)
            return self

        return method

    def run(self):
        kmer_graph = self.serializer_driver.run()
        the_serializer = (serializer
                          .Serializer(kmer_graph, colors=self.serializer_driver.retriever.colors))
        the_serializer.to_json_serializable()
        return CollapsedKmerUnitgGraphExpectation(the_serializer.unitig_graph)


class TestCollapseKmerUnitigsCreatesSingleUnitig(object):
    def test_with_missing_kmer(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_contig_to_retrieve('GTT'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('GTT').has_coverages(0, 1)

    def test_with_one_kmer(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAC', 1)
                  .with_contig_to_retrieve('GTT'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('GTT').has_coverages([1, 1])

    def test_with_two_linked_kmers(self):
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.....C..')
                  .with_kmer('AAC', 0, 'a.......')
                  .with_contig_to_retrieve('AAAC'))

        # when
        expect = driver.run()

        # then
        assert expect.has_kmers('AAAC')

    def test_with_three_linked_kmers(self):
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.....C..')
                  .with_kmer('AAC', 0, 'a....C..')
                  .with_kmer('ACC', 0, 'a.......')
                  .with_contig_to_retrieve('AAACC'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('AAACC')


class TestCollapseKmerUnitigs(object):
    def test_with_two_unlinked_missing_kmers_creates_single_unitig(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA')
                  .with_kmer('AAC')
                  .with_contig_to_retrieve('GTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GTTG').has_coverages([0, 1], [0, 1])
        expect.has_n_nodes(1)
        expect.has_n_edges(0)

    def test_with_two_unlinked_kmers_creates_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 1)
                  .with_kmer('AAC', 1)
                  .with_contig_to_retrieve('GTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_one_node_with_repr('GTT').has_coverages(1, 1)
        expect.has_one_node_with_repr('G').has_coverages(1, 1)
        expect.has_n_edges(1)

    def test_with_two_node_path_and_three_node_cycle_results_in_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.....G.')
                  .with_kmer('ACG', 1, 'a.g.A...')
                  .with_kmer('CGA', 1, 'a....C..')
                  .with_kmer('GAC', 1, '.c....G.')
                  .with_contig_to_retrieve('AAACGAC'))

        # when
        expect = driver.run()

        # then
        expect.has_kmers('AAAC', 'GAC')
        expect.has_one_node_with_repr('AAAC').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('GAC').has_coverages([1, 1], [1, 1], [1, 1])
        expect.has_n_edges(3)
        expect.has_n_missing_edges(0)

    def test_four_node_path_with_one_node_bubble_in_three_nodes(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAC', 1, '.....C..')
                  .with_kmer('ACC', 1, 'a....CG.')
                  .with_kmer('CCC', 1, 'a.....G.')
                  .with_kmer('CCG', 1, 'ac..A...')
                  .with_kmer('CGA', 1, '.c......')
                  .with_contig_to_retrieve('AACCCCGA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('AACC').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('C').has_coverages([1, 1])
        expect.has_one_node_with_repr('GA').has_coverages([1, 1])
        expect.has_n_nodes(3)

    def test_two_kmers_one_kmer_apart_do_not_collapse(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('GAA', 1, '........')
                  .with_kmer('ACC', 1, '........')
                  .with_contig_to_retrieve('GGTTC'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(3)
        expect.has_n_edges(2)
        expect.has_one_node_with_repr('GGT').has_coverages([1, 1])
        expect.has_one_node_with_repr('T').has_coverages([0, 1])
        expect.has_one_node_with_repr('C').has_coverages([1, 1])

    def test_two_linked_kmers_with_incoming_edge_and_missing_kmer_returns_three_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.c...C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .with_contig_to_retrieve('GTTTA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('A').has_coverages([0, 1])
        expect.has_one_node_with_repr('G').has_coverages([0, 0])
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])

        expect.has_n_nodes(3)
        expect.has_n_edges(2)


class TestToJson(object):
    def test_two_linked_kmers_are_jsonifiable(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2)
                         .with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
                         .with_kmer('AAC', [1, 0], ['a.......', '........']))
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        kmer_json = serializer.Serializer(kmer_graph, colors=[0, 1, 2]).to_json()

        # then
        json.loads(kmer_json)  # does not raise


class TestToJsonSerializable(object):
    def test_two_linked_kmers_collapse_to_one_kmer(self):
        # given
        driver = (JsonSerializableTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .with_contig_to_retrieve('GTTT'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1).has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])

    def test_two_kmers_one_kmer_apart_do_not_collapse(self):
        # given
        driver = (JsonSerializableTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 1, '........')
                  .with_kmer('ACC', 1, '........')
                  .with_contig_to_retrieve('GGTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GGT').has_coverages([1, 1])
        expect.has_one_node_with_repr('T').has_coverages([0, 1])
        expect.has_one_node_with_repr('G').has_coverages([1, 1])

        expect.has_n_nodes(3)
        expect.has_n_edges(2)

    def test_unlinked_kmers_followed_by_two_linked_kmers_collapse_to_two_unitigs(self):
        # given
        driver = (JsonSerializableTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .with_contig_to_retrieve('GTTTAA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('AA').has_coverages([0, 1], [0, 1])
        expect.has_n_nodes(2)
        expect.has_n_edges(1)

    def test_linked_kmers_with_outgoing_edge_surrounded_by_missing_kmers(self):
        # given
        driver = (JsonSerializableTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a....C..')
                  .with_contig_to_retrieve('CAAGTTTAGG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('TT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('CAAGT').has_coverages([[0, 1] for _ in range(3)])
        expect.has_one_node_with_repr('GGT').has_coverages([0, 0])
        expect.has_one_node_with_repr('AGG').has_coverages([0, 1])
        expect.has_n_nodes(4)
        expect.has_n_edges(3)
