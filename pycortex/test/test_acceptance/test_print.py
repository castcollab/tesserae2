import contextlib
import io
from subprocess import check_call, check_output

import pytest

from pycortex.__main__ import main
from pycortex.test.builders.mccortex_builder import MccortexFactory, MCCORTEX

PYCORTEX_COMMAND = ['python', '-m', 'pycortex']


class TestPrintCommand(object):
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

        assert expected_kmer == pycortex_output.getvalue().rstrip()

    @pytest.mark.xfail
    def test_prints_three_kmers(self, tmpdir):
        # given
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', 'ACCAA')
                   .with_kmer_size(3))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])

        expected_kmers = [
            'ACC: ACC 1 ....A...',
            'CAA: CAA 1 .c......',
            'CCA: CCA 1 a...A...',
        ]

        # when

        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['print', '--graph', output_graph])

        assert expected_kmers == pycortex_output.getvalue().split()
