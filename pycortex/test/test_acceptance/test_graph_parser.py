import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import check_call

from os.path import join

from pycortex.cortex_graph import CortexGraphHeader

BIN_DIR = os.environ['BIN_DIR']
MCCORTEX = join(BIN_DIR, 'mccortex31')


class TestCortexGraphParsing(object):
    def test_parses_a_graph_header(self, tmpdir):
        sample_name = b'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3
        input_fasta = tmpdir.join('input.fasta')
        output_graph = tmpdir.join('output.ctx')
        input_fasta.write(SeqRecord(Seq(dna_sequence)).format('fasta'))

        check_call(
            [MCCORTEX, 'build',
             '--kmer', str(kmer_size),
             '--sample', sample_name,
             '-1', str(input_fasta),
             output_graph])

        with open(output_graph, 'rb') as fh:
            cgh = CortexGraphHeader.from_stream(fh)

        expected = CortexGraphHeader(version=6,
                                     kmer_size=kmer_size,
                                     kmer_container_size=1,
                                     num_colors=1,
                                     mean_read_lengths=(len(dna_sequence),),
                                     mean_total_sequence=(len(dna_sequence),),
                                     sample_names=(sample_name,))
        assert cgh == expected
