import cortexpy.test.driver.command as command


class TestSubgraphs(object):
    def test_traverses_single_subgraph_as_single_subgraph(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_n_graphs(1)

    def test_traverses_two_subgraphs_as_two_subgraphs(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC', 'AAA')
        d.with_subgraph_output()
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_n_graphs(2)
