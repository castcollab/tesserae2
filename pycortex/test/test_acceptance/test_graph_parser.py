from subprocess import check_call

from pycortex.graph.parser.header import Header, header_from_stream
from pycortex.graph.parser.random_access import RandomAccess
from pycortex.graph.parser.streaming import kmer_generator_from_stream
from pycortex.test.builders.graph.body import KmerRecord, as_edge_set, print_kmer
from pycortex.test.builders.mccortex_builder import MCCORTEX, MccortexFactory


class TestCortexGraph(object):
    def test_parses_a_graph_header(self, tmpdir):
        # given
        sample_name = b'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3

        factory = (MccortexFactory()
                   .with_dna_sequence(sample_name, dna_sequence)
                   .with_kmer_size(kmer_size))

        expected_header = Header(version=6,
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

        header = header_from_stream(open(output_graph, 'rb'))

        # then
        assert header == expected_header

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
        kmer_generator = kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
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
        cg = RandomAccess(open(output_graph, 'rb'))

        # when
        actual = cg.get_kmer('AAC')

        # then
        print_kmer(actual)

        assert actual.kmer == expected.kmer
        assert actual.coverage == expected.coverage
        assert actual.edges == expected.edges

    def test_parses_a_graph_with_kmer_size_9(self, tmpdir):
        # given
        factory = (MccortexFactory()
                   .with_dna_sequence(b'sample_0', 'ACGTTCCCC')
                   .with_kmer_size(9))

        expected_kmers = [
            KmerRecord('ACGTTCCCC', (1,), (as_edge_set('........'),)),
        ]

        # when
        output_graph = factory.build(tmpdir)

        check_call([MCCORTEX, 'view', '-k', output_graph])

        kmer_generator = kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            print_kmer(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges
