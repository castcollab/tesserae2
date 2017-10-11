import contextlib
import io
import json

import attr
from pycortex.__main__ import main
import pycortex.test.builder as builder
import pycortex.test.runner as runner


@attr.s(slots=True)
class PycortexPrintOutputParser(object):
    output = attr.ib()

    def get_kmer_strings(self):
        return self.output.rstrip().split('\n')


class Test(object):
    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'AAG 1 ......G.',
            'ACC 1 .......T',
            'AGG 1 a......T',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()


class TestTermWithRecord(object):
    def test_prints_single_kmer(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('ACCAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CAA: CAA 1 .c......'

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', 'CAA'])

        # then
        assert [expected_kmer] == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_one_missing_missing_kmer(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('AAAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CCC: GGG missing'

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', 'GGG'])

        # then
        assert [expected_kmer] == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers(self, tmpdir):
        # given
        record = 'ACCAA'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size).build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 ....A...',
            'CCA: CCA 1 a...A...',
            'CAA: CAA 1 .c......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 .......T',
            'AGG: CCT 1 A......t',
            'AAG: CTT 1 .C......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_one_missing_and_one_revcomp_kmer(self, tmpdir):
        # given
        record = 'ACCTT'
        search_record = 'ACTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACT: ACT missing',
            'AAG: CTT 1 .C......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', search_record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()


class TestOutputTypeJSON(object):
    def test_outputs_json(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        # when
        completed_process = (runner
                             .Pycortex(True)
                             .view(['--record', record, '--output-type', 'json', output_graph]))
        stdout = completed_process.stdout.decode()

        # then
        assert completed_process.returncode == 0, completed_process
        graph = json.loads(stdout)
        assert graph['directed']
