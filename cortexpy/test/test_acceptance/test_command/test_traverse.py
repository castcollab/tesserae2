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

    def test_traverses_all_colors_without_color_specification(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_record('CCCG', name='0')
        d.with_record('CCGAA', name='1')
        d.with_initial_contigs('CCC')
        d.without_traversal_colors()
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_nodes('CCC', 'CCG', 'CGA', 'GAA')


class TestLogging(object):
    def test_without_logging_args_emits_info_log_message(self, tmpdir, caplog):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)

        # when
        d.run()

        # then
        assert 'Log level is' in caplog.text

    def test_with_logging_args_emits_info_log_message(self, tmpdir, caplog):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)
        d.with_verbose_arg()
        d.with_logging_interval(0)

        # when
        d.run()

        # then
        assert 'Log level is' in caplog.text
        assert 'current graph size: 3' in caplog.text

    def test_with_logging_interval_2_does_not_report_current_graph_size(self, tmpdir, caplog):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)
        d.with_verbose_arg()
        d.with_logging_interval(10)

        # when
        d.run()

        # then
        assert 'Log level is' in caplog.text
        assert 'current graph size: 3' not in caplog.text
