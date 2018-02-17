import networkx as nx

from cortexpy.graph import interactor


def test_only_follows_one_color_with_color_specified():
    # given
    graph = nx.MultiDiGraph()
    graph.add_edge('AAA', 'AAT', key=0)
    graph.add_edge('AAT', 'ATA', key=1)

    # when
    paths = list(interactor.Contigs(graph, color=0).all_simple_paths())

    # then
    assert ['AAAT'] == [str(p.seq) for p in paths]


def test_follows_two_colors_with_no_color_specified():
    # given
    graph = nx.MultiDiGraph()
    graph.add_edge('AAA', 'AAT', key=0)
    graph.add_edge('AAT', 'ATA', key=1)

    # when
    paths = list(interactor.Contigs(graph).all_simple_paths())

    # then
    assert {'AAATA'} == set([str(p.seq) for p in paths])


def test_follows_three_colors_with_no_color_specified():
    # given
    graph = nx.MultiDiGraph()
    graph.add_edge('AAA', 'AAT', key=0)
    graph.add_edge('AAT', 'ATA', key=1)
    graph.add_edge('AAT', 'ATC', key=2)

    # when
    paths = list(interactor.Contigs(graph).all_simple_paths())

    # then
    assert {'AAATA', 'AAATC'} == set([str(p.seq) for p in paths])
