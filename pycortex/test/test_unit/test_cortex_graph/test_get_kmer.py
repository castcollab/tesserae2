from io import BytesIO

from hypothesis import given
from hypothesis import strategies as s
from math import ceil

from pycortex.cortex_graph import CortexGraphRandomAccessParser
from pycortex.test.builders.graph_body_builder import kmers
from pycortex.test.builders.graph_builder import CortexGraphBuilder


class TestGetKmer(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=0, max_value=5))
    def test_get_record(self, data, kmer_size, num_colors, n_kmers):
        # given
        builder = CortexGraphBuilder()
        builder.with_kmer_size(kmer_size)
        builder.with_num_colors(num_colors)
        builder.sort_kmers = True

        expected_kmers = []
        seen_kmers = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            while kmer.kmer in seen_kmers:
                kmer = data.draw(kmers(kmer_size, num_colors))
            seen_kmers.add(kmer.kmer)
            builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        cg = CortexGraphRandomAccessParser(builder.build())

        # when
        for expected_kmer in expected_kmers:
            kmer = cg.get_kmer(expected_kmer.kmer)
            assert kmer

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert expected_kmer.coverage == kmer.coverage
            assert expected_kmer.edges == kmer.edges
