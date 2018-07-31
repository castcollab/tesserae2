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
        expect.has_n_nodes(3)

    def test_traverses_two_subgraphs_as_single_graph(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC', 'TTT')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_nodes('CCC', 'CCG', 'CGC', 'AAA')

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
    def test_without_logging_args_emits_info_log_message(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)

        # when
        stderr = d.run_for_stderr()

        # then
        assert 'Log level is INFO' in stderr

    def test_with_logging_args_emits_info_log_message(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)
        d.with_verbose_arg()
        d.with_logging_interval(0)

        # when
        stderr = d.run_for_stderr()

        # then
        assert 'Log level is DEBUG' in stderr
        assert 'current graph size: 3' in stderr

    def test_with_logging_interval_2_does_not_report_current_graph_size(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)
        d.with_verbose_arg()
        d.with_logging_interval(10)

        # when
        stderr = d.run_for_stderr()

        # then
        assert 'Log level is DEBUG' in stderr
        assert 'current graph size: 3' not in stderr

    def test_with_silent_arg_emits_no_info_log_message(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC')
        d.with_kmer_size(3)
        d.with_silent_arg()

        # when
        stderr = d.run_for_stderr()

        # then
        assert 'Log level is' not in stderr


class Test(object):
    def test_creates_two_transcripts_from_four_records_in_four_colors_in_four_graphs(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        nodes = ['CCCG', 'CCGC', 'AAAT', 'AATA']
        for idx, seq in enumerate(nodes):
            if idx != 0:
                d.with_extra_graph()
            d.with_record(seq)

        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_nodes('CCC', 'CCG', 'CGC', 'AAA', 'AAT', 'ATA')
