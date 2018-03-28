from hypothesis import given
from hypothesis import strategies as s

from cortexpy.graph.parser.streaming import kmer_list_generator_from_stream_and_header
from cortexpy.test.builder.graph.body import Body, KmerRecord
from cortexpy.test.builder.graph.kmer import kmers
from cortexpy.test.mock.graph import Header


class TestStreamKmerGenerator(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=260),
           s.integers(min_value=0, max_value=10),
           s.integers(min_value=0, max_value=4))
    def test_parses_records(self, data, kmer_size, num_colors, n_kmers):
        # given
        builder = Body(kmer_size=kmer_size)

        expected_kmers = []
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        header = Header(kmer_size, builder.kmer_container_size, num_colors)

        # when
        for kmer_list, expected_kmer in zip(
            kmer_list_generator_from_stream_and_header(builder.build(), header),
            expected_kmers
        ):
            # then
            assert expected_kmer.kmer == ''.join(kmer_list)

    def test_parses_aac_kmer(self):
        kmer_container_size = 1
        kmer_size = 3
        num_colors = 0
        header = Header(kmer_size, kmer_container_size, num_colors)
        builder = Body(kmer_container_size)

        expected_kmer = KmerRecord('AAC', [], [])
        builder.with_kmer_record(expected_kmer)

        kmer_list = next(kmer_list_generator_from_stream_and_header(builder.build(), header))

        assert expected_kmer.kmer == ''.join(kmer_list)
