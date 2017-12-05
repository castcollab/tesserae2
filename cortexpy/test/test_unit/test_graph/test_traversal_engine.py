import attr
import pytest

import cortexpy as cp
import cortexpy.test.builder as builder
import cortexpy.test.expectation as expectation


@attr.s(slots=True)
class TraversalEngineTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def run(self):
        graph = cp.graph.TraversalEngine(self.graph_builder.build()).traverse()
        return expectation.kmer.GraphExpectation(graph)


class Test(object):
    @pytest.mark.skip('Not implemented yet')
    def test_on_empty_graph_returns_empty_graph(self):
        # given
        driver = TraversalEngineTestDriver().with_kmer_size(3)

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_nodes(0)
         .has_n_edges(0))
