import os
from itertools import repeat

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import check_call, check_output

from os.path import join

from attr import Factory

from pycortex.cortex_graph import CortexGraphHeader, CortexGraphStreamingParser, \
    CortexGraphRandomAccessParser
from pycortex.test.builders.graph_body_builder import KmerRecord

BIN_DIR = os.environ['BIN_DIR']
MCCORTEX = join(BIN_DIR, 'mccortex31')


def as_edge_set(edge_set_string):
    assert len(edge_set_string) == 8
    return tuple([edge != '.' for edge in edge_set_string])


@attr.s
class MccortexFactory(object):
    kmer_size = attr.ib(3)
    sequences = attr.ib(Factory(list))

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_dna_sequence(self, name, sequence):
        self.sequences.append((name, sequence))
        return self

    def build(self, tmpdir):
        command = [MCCORTEX, 'build', '--sort', '--kmer', str(self.kmer_size)]
        input_fasta = str(tmpdir.join('input.fasta'))
        with open(input_fasta, 'w') as fh:
            for name, dna_sequence in self.sequences:
                fh.write(SeqRecord(Seq(dna_sequence)).format('fasta'))
                command.extend(['--sample', name, '-1', input_fasta])

        output_graph = str(tmpdir.join('output.ctx'))
        command.append(output_graph)

        check_call(command)

        return output_graph


def print_kmer(kmer):
    print(kmer)
    print(kmer.kmer)
    print(kmer.coverage)
    print(kmer.edges)


class TestCortexGraph(object):
    def test_parses_a_graph_header(self, tmpdir):
        # given
        sample_name = b'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3

        factory = (MccortexFactory()
                   .with_dna_sequence(sample_name, dna_sequence)
                   .with_kmer_size(kmer_size))

        expected_header = CortexGraphHeader(version=6,
                                            kmer_size=kmer_size,
                                            record_size=13,
                                            kmer_container_size=1,
                                            num_colors=1,
                                            mean_read_lengths=(len(dna_sequence),),
                                            mean_total_sequence=(len(dna_sequence),),
                                            sample_names=(sample_name,))

        # when
        output_graph = factory.build(tmpdir)

        check_call([MCCORTEX, 'view', '-k', output_graph])

        cg = CortexGraphStreamingParser(open(output_graph, 'rb'))

        # then
        assert cg.header == expected_header

    def test_parses_a_graph(self, tmpdir):
        # given
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', 'ACGTT')
                   .with_kmer_size(3))

        expected_kmers = [
            KmerRecord('AAC', (1,), (as_edge_set('......G.'),)),
            KmerRecord('ACG', (2,), (as_edge_set('A......T'),)),
        ]

        # when
        output_graph = factory.build(tmpdir)

        check_call([MCCORTEX, 'view', '-k', output_graph])

        cg = CortexGraphStreamingParser(open(output_graph, 'rb'))

        # then
        actual_kmers = list(cg.kmers())
        for kmer in actual_kmers:
            print_kmer(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges

    def test_retrieves_kmer_by_random_access(self, tmpdir):
        # given
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', 'ACGTTT')
                   .with_kmer_size(3))

        expected = KmerRecord('AAC', (1,), (as_edge_set('A.....G.'),))
        output_graph = factory.build(tmpdir)
        check_call([MCCORTEX, 'view', '-k', output_graph])
        cg = CortexGraphRandomAccessParser(open(output_graph, 'rb'))

        # when
        actual = cg.get_kmer('AAC')

        # then
        print_kmer(actual)

        assert actual.kmer == expected.kmer
        assert actual.coverage == expected.coverage
        assert actual.edges == expected.edges
