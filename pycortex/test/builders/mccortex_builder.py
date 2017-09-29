import os
from os.path import join
from subprocess import check_call

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory

BIN_DIR = os.environ['BIN_DIR']
MCCORTEX = join(BIN_DIR, 'mccortex31')


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
