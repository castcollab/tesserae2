import pytest
from hypothesis import given
from hypothesis import strategies as s

from pycortex.cortex_graph import CortexGraphRandomAccessParser, CortexGraphRandomAccessError
from pycortex.test.builders.graph_body_builder import kmers, KmerRecord, as_edge_set
from pycortex.test.builders.graph_builder import CortexGraphBuilder


class TestGetKmer(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=0, max_value=5))
    def test_get_records(self, data, kmer_size, num_colors, n_kmers):
        # given
        builder = CortexGraphBuilder()
        builder.with_kmer_size(kmer_size)
        builder.with_num_colors(num_colors)
        builder.with_sorted_kmers()

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

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert expected_kmer.coverage == kmer.coverage
            assert expected_kmer.edges == kmer.edges

    def test_raises_on_missing_kmer(self):
        # given
        builder = CortexGraphBuilder()
        builder.with_kmer_size(3)
        builder.with_sorted_kmers()

        cg = CortexGraphRandomAccessParser(builder.build())

        # when
        with pytest.raises(CortexGraphRandomAccessError):
            cg.get_kmer('AAA')


class TestGetKmerForString(object):
    def test_gets_AAA_for_TTT_query(self):
        # given
        builder = CortexGraphBuilder()
        builder.with_kmer_size(3)
        builder.with_num_colors(1)
        builder.with_sorted_kmers()

        expected_kmer = KmerRecord('AAA', [1], [as_edge_set('........')])
        builder.with_kmer_record(expected_kmer)

        cg = CortexGraphRandomAccessParser(builder.build())

        # when
        assert expected_kmer.kmer == cg.get_kmer_for_string('AAA').kmer
        assert expected_kmer.kmer == cg.get_kmer_for_string('TTT').kmer
