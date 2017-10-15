import attr

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
        assert expect.has_kmers('AAACC')


@attr.s
class KmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages(self, *coverages):
        assert self.kmer_node['coverage'] == list(coverages)
        return self

    def is_missing(self):
        assert self.kmer_node['is_missing']
        return self


@attr.s(slots=True)
class CollapsedKmerUnitgGraphExpectation(object):
    graph = attr.ib()
    nodes_by_repr = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.nodes_by_repr = {data['repr']: node for node, data in self.graph.nodes.data()}

    def has_n_kmers(self, n):
        assert len(self.graph) == n
        return self

    def has_kmers(self, *kmer_reprs):
        assert set(self.nodes_by_repr.keys()) == set(kmer_reprs)
        return self

    def has_kmer(self, kmer_repr):
        assert kmer_repr in self.nodes_by_repr
        return KmerNodeExpectation(self.graph.node[self.nodes_by_repr[kmer_repr]])


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
        expect.has_kmers('AACC', 'C', 'GA')


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
        expect.has_kmer('TTT').has_coverages(1, 1)
        expect.has_kmer('GTT').has_coverages(1, 0)
        expect.has_kmer('ATT').is_missing()


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
        assert kmer_json.startswith('{"directed')


class TestToJsonSerializable(object):
    def test_two_linked_kmers_have_node_id_0(self):
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
        expect.has_n_kmers(1).has_kmer('GTTT')
