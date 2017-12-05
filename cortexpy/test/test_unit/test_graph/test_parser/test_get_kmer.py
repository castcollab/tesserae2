import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as s
import cortexpy.test.builder as builder
import cortexpy.graph.parser as parser
from cortexpy.test.builder.graph.body import KmerRecord, as_edge_set
from cortexpy.test.builder.graph.kmer import kmers


class TestDunderGetitemDunder(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=0, max_value=5))
    def test_record_retrieval(self, data, kmer_size, num_colors, n_kmers):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(kmer_size)
                         .with_num_colors(num_colors))

        expected_kmers = []
        seen_kmers = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            while kmer.kmer in seen_kmers:
                kmer = data.draw(kmers(kmer_size, num_colors))
            seen_kmers.add(kmer.kmer)
            graph_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        cg = parser.RandomAccess(graph_builder.build())

        # when
        for expected_kmer in expected_kmers:
            kmer = cg[expected_kmer.kmer]

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    def test_raises_on_missing_kmer(self):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(3)

        cg = parser.RandomAccess(graph_builder.build())

        # when
        with pytest.raises(KeyError):
            cg['AAA']


class TestGetKmerForString(object):
    def test_gets_aaa_for_ttt_query(self):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(3)
        graph_builder.with_num_colors(1)

        expected_kmer = KmerRecord('AAA', [1], [as_edge_set('........')])
        graph_builder.with_kmer_record(expected_kmer)

        cg = parser.RandomAccess(graph_builder.build())

        # when
        assert expected_kmer.kmer == cg.get_kmer_for_string('AAA').kmer
        assert expected_kmer.kmer == cg.get_kmer_for_string('TTT').kmer


class TestDunderIterDunder(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=10),
           s.integers(min_value=1, max_value=10),
           s.integers(min_value=0, max_value=4))
    def test_parses_records(self, data, kmer_size, num_colors, n_kmers):
        # given
        graph_builder = (builder.Graph()
                         .without_sorted_kmers()
                         .with_kmer_size(kmer_size)
                         .with_num_colors(num_colors))

        expected_kmers = []
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            graph_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)
        ra_parser = parser.RandomAccess(graph_builder.build())

        # when
        for kmer, expected_kmer in zip(ra_parser, expected_kmers):
            # then
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    def test_gets_aaa(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(1))

        expected_kmer = KmerRecord('AAA', (1,), [as_edge_set('........')])
        graph_builder.with_kmer_record(expected_kmer)

        cg = parser.RandomAccess(graph_builder.build())

        # when
        for kmer in cg:
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    def test_works_with_no_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(1))

        cg = parser.RandomAccess(graph_builder.build())

        # when/then
        assert list(cg) == []
