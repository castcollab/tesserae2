import os
from itertools import repeat

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import check_call, check_output

from os.path import join

from pycortex.cortex_graph import CortexGraphHeader, CortexGraphFromStream
from pycortex.test.builders.graph_body_builder import KmerRecord

BIN_DIR = os.environ['BIN_DIR']
MCCORTEX = join(BIN_DIR, 'mccortex31')


def as_edge_set(edge_set_string):
    assert len(edge_set_string) == 8
    return tuple([edge != '.' for edge in edge_set_string])


class TestCortexGraphParsing(object):
    def test_parses_a_graph_header(self, tmpdir):
        # given
        sample_name = b'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3
        input_fasta = tmpdir.join('input.fasta')
        output_graph = str(tmpdir.join('output.ctx'))
        input_fasta.write(SeqRecord(Seq(dna_sequence)).format('fasta'))

        expected_header = CortexGraphHeader(version=6,
                                            kmer_size=kmer_size,
                                            kmer_container_size=1,
                                            num_colors=1,
                                            mean_read_lengths=(len(dna_sequence),),
                                            mean_total_sequence=(len(dna_sequence),),
                                            sample_names=(sample_name,))

        expected_kmers = [
            KmerRecord(tuple('AAC'), (1,), (as_edge_set('......G.'),)),
            KmerRecord(tuple('ACG'), (2,), (as_edge_set('A......T'),)),
        ]

        # when
        check_call(
            [MCCORTEX, 'build',
             '--sort',
             '--kmer', str(kmer_size),
             '--sample', sample_name,
             '-1', str(input_fasta),
             output_graph])

        check_call([MCCORTEX, 'view', '-k', output_graph])

        cg = CortexGraphFromStream(open(output_graph, 'rb'))

        # then
        assert cg.header == expected_header
        actual_kmers = list(cg.kmers())
        for kmer in actual_kmers:
            print(kmer)
            print(kmer.kmer)
            print(kmer.coverage)
            print(kmer.edges)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges
