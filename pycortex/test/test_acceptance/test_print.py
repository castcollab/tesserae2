import contextlib
import io
from subprocess import check_call

import attr
from pycortex.__main__ import main
from pycortex.test.builders.mccortex_builder import MccortexFactory, MCCORTEX

PYCORTEX_COMMAND = ['python', '-m', 'pycortex']


@attr.s(slots=True)
class PycortexPrintOutputParser(object):
    output = attr.ib()

    def get_kmer_strings(self):
        return self.output.rstrip().split('\n')


class TestPrintCommandWithRecord(object):
    def test_prints_single_kmer(self, tmpdir):
        # given
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', 'ACCAA')
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmer = 'CAA 1 .c......'

        # when

        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph, '--record', 'CAA'])

        assert [expected_kmer] == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers(self, tmpdir):
        # given
        record = 'ACCAA'
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', record)
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmers = [
            'ACC: ACC 1 ....A...',
            'CCA: CCA 1 a...A...',
            'CAA: CAA 1 .c......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', record)
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmers = [
            'ACC: ACC 1 .......T',
            'AGG: CCT 1 A......t',
            'AAG: CTT 1 .C......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_one_missing_and_one_revcomp_kmer(self, tmpdir):
        # given
        record = 'ACCTT'
        search_record = 'ACTT'
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', record)
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmers = [
            'ACT: ACT missing',
            'AAG: CTT 1 .C......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph, '--record', search_record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()


class TestPrintCommand(object):
    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', record)
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmers = [
            'AAG 1 ......G.',
            'ACC 1 .......T',
            'AGG 1 a......T',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()
