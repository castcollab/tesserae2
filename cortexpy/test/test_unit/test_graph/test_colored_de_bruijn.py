from hypothesis import given, strategies as s

from cortexpy.graph.interactor import make_multi_graph, Interactor
from cortexpy.test.builder.graph.colored_de_bruijn import get_cdb_builder


class Test(object):
    def test_with_2_nodes_has_len_2(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 .....C..')
        b.with_kmer('AAC 0 a.......')

        # when
        cdb = b.build()
        ucdb = make_multi_graph(cdb)

        # then
        assert 2 == len(cdb)
        assert 2 == len(ucdb)


class TestNodes(object):
    def test_nodes_returns_two_nodes(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 .....C..')
        b.with_kmer('AAC 0 a.......')

        # when
        cdb = b.build()

        # then
        assert {'AAA', 'AAC'} == set(cdb.nodes())

    @given(s.booleans())
    def test_remove_node(self, remove_first):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 .....C..')
        b.with_kmer('AAC 0 a.......')

        cdb = b.build()

        # when
        if remove_first:
            cdb.remove_node('AAC')

            # then
            assert {'AAA'} == set(cdb.nodes())
        else:
            cdb.remove_node('AAA')

            # then
            assert {'AAC'} == set(cdb.nodes())

    @given(s.lists(s.sampled_from(('AAC', 'AAA', 'ACC'))))
    def test_remove_nodes(self, removal_nodes):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 .....C..')
        b.with_kmer('AAC 0 a....C..')
        b.with_kmer('ACC 0 a.......')
        nodes = {'AAC', 'AAA', 'ACC'}

        cdb = b.build()

        # when
        cdb.remove_nodes_from(removal_nodes)

        # then
        assert (nodes - set(removal_nodes)) == set(cdb.nodes())


class TestEdges(object):
    def test_neighbors_of_three_colors_returns_colors_0_and_2(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 0 0 .....C.. ........ .....C..')
        b.with_kmer('AAC 0 0 0 a....... ........ a.......')

        # when
        cdb = b.build()

        # then
        assert {0, 2} == set(cdb.succ['AAA']['AAC'])
        assert set() == set(cdb.pred['AAA']['AAC'])
        assert {0, 2} == set(cdb['AAA']['AAC'])

        assert set() == set(cdb.succ['AAC']['AAA'])
        assert {0, 2} == set(cdb.pred['AAC']['AAA'])
        assert set() == set(cdb['AAC']['AAA'])


class TestInOutEdges(object):
    def test_single_kmer(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 ......G.')
        b.with_kmer('AAG 0 a.......')
        cdb = b.build()
        seed = 'AAA'

        # when / then
        assert [('AAA', 'AAG')] == list(cdb.out_edges(seed))
        assert [] == list(cdb.in_edges(seed))

    def test_single_kmer_revcomp(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 ......G.')
        b.with_kmer('AAG 0 a.......')
        cdb = b.build()
        seed = 'TTT'

        # when / then
        assert [('CTT', 'TTT')] == list(cdb.in_edges(seed))
        assert [] == list(cdb.out_edges(seed))


class TestConsistentCortexGraph(object):
    @given(s.sampled_from(('AAA', 'TTT')))
    def test_single_kmer_revcomp_seed(self, seed):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAA 0 ......G.')
        b.with_kmer('AAG 0 a.......')
        cdb = b.build()

        # when
        graph = Interactor(cdb, colors=None).make_graph_nodes_consistent([seed]).graph

        # then
        if seed == 'AAA':
            assert [] == list(graph.in_edges(seed))
            assert [('AAA', 'AAG')] == list(graph.out_edges(seed))
        else:
            assert [('CTT', 'TTT')] == list(graph.in_edges(seed))
            assert [] == list(graph.out_edges(seed))

    def test_gets_correct_neighbors_of_kmer(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('AAC 0 .......T')
        b.with_kmer('ACT 0 a.....G.')
        b.with_kmer('CAG 0 .......T')
        cdb = b.build()
        seed = 'AAC'

        # when
        graph = Interactor(cdb, colors=None).make_graph_nodes_consistent([seed]).graph

        # then
        assert ['CTG'] == list(graph['ACT'])
        assert ['CTG'] == list(graph.succ['ACT'])
        assert ['AAC'] == list(graph.pred['ACT'])
        assert [('ACT', 'CTG')] == list(graph.out_edges('ACT'))
        assert [('AAC', 'ACT')] == list(graph.in_edges('ACT'))

        assert [] == list(graph['CTG'])
        assert [] == list(graph.succ['CTG'])
        assert ['ACT'] == list(graph.pred['CTG'])
        assert [] == list(graph.out_edges('CTG'))
        assert [('ACT', 'CTG')] == list(graph.in_edges('CTG'))
