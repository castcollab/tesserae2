import json

import attr

import cortexpy.graph as graph
import cortexpy.graph.serializer as serializer
import cortexpy.test.builder as builder
from cortexpy.graph import traversal, parser
from cortexpy.test.expectation.kmer import CollapsedKmerUnitgGraphExpectation


@attr.s(slots=True)
class SerializerTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    contig_to_retrieve = attr.ib(None)
    retriever = attr.ib(None)
    traverse = attr.ib(False)
    retrieve = attr.ib(False)
    traversal_start_kmer = attr.ib(None)
    traversal_colors = attr.ib((0,))

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def traverse_with_start_kmer_and_color(self, start_kmer, color):
        self.traverse = True
        self.traversal_start_kmer = start_kmer
        self.traversal_colors = [color]
        return self

    def retrieve_contig(self, contig):
        self.retrieve = True
        self.contig_to_retrieve = contig
        return self

    def run(self):
        if self.retrieve:
            self.retriever = graph.ContigRetriever(self.graph_builder.build())
            return self.retriever.get_kmer_graph(self.contig_to_retrieve)
        elif self.traverse:
            traverser = traversal.Engine(parser.RandomAccess(self.graph_builder.build()),
                                         traversal_colors=self.traversal_colors)
            return traverser.traverse_from(self.traversal_start_kmer).graph
        else:
            raise Exception("Need to load a command")


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
                     .UnitigCollapser(kmer_graph)
                     .collapse_kmer_unitigs())
        return CollapsedKmerUnitgGraphExpectation(collapser.unitig_graph)


class TestCollapseKmerUnitigsCreatesSingleUnitig(object):
    def test_with_missing_kmer(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .retrieve_contig('GTT'))

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
                  .retrieve_contig('GTT'))

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
                  .retrieve_contig('AAAC'))

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
                  .retrieve_contig('AAACC'))

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
                  .retrieve_contig('GTTG'))

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
                  .retrieve_contig('GTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_one_node_with_repr('GTT').has_coverages(1, 1)
        expect.has_one_node_with_repr('G').has_coverages(1, 1)
        expect.has_n_edges(1)

    def test_two_linked_kmers_collapse_to_one_kmer(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .retrieve_contig('GTTT'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1).has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])

    def test_with_two_node_path_and_three_node_cycle_results_in_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.....G.')
                  .with_kmer('ACG', 1, 'a.g.A...')
                  .with_kmer('CGA', 1, 'a....C..')
                  .with_kmer('GAC', 1, '.c....G.')
                  .retrieve_contig('AAACGAC'))

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
                  .retrieve_contig('AACCCCGA'))

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
                  .retrieve_contig('GGTTC'))

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
                  .retrieve_contig('GTTTA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('A').has_coverages([0, 1])
        expect.has_one_node_with_repr('G').has_coverages([0, 0])
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])

        expect.has_n_nodes(3)
        expect.has_n_edges(2)

    def test_unlinked_kmers_followed_by_two_linked_kmers_collapse_to_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .retrieve_contig('GTTTAA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('AA').has_coverages([0, 1], [0, 1])
        expect.has_n_nodes(2)
        expect.has_n_edges(1)

    def test_linked_kmers_with_outgoing_edge_surrounded_by_missing_kmers_returns_four_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a....C..')
                  .retrieve_contig('CAAGTTTAGG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('TT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('CAAGT').has_coverages([[0, 1] for _ in range(3)])
        expect.has_one_node_with_repr('GGT').has_coverages([0, 0])
        expect.has_one_node_with_repr('AGG').has_coverages([0, 1])
        expect.has_n_nodes(4)
        expect.has_n_edges(3)


class TestToJson(object):
    def test_two_linked_kmers_are_jsonifiable(self):
        # given
        color_names = 'samp1', 'samp2'
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2)
                         .with_color_names(*color_names)
                         .with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
                         .with_kmer('AAC', [1, 0], ['a.......', '........']))
        retriever = graph.ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        kmer_json = serializer.Serializer(kmer_graph).to_json()

        # then
        kmer_data = json.loads(kmer_json)  # does not raise
        assert kmer_data['graph']['colors'] == list(range(3))
        assert kmer_data['graph']['sample_names'] == list(color_names) + ['retrieved_contig']
