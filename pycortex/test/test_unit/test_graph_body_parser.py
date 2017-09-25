from io import BytesIO

import attr
from hypothesis import given, example
from hypothesis import strategies as s
from math import ceil

from pycortex.cortex_graph import kmer_generator_from_stream
from pycortex.test.builders.graph_body_builder import CortexGraphBodyBuilder, KmerRecord


@s.composite
def kmers(draw, kmer_size, num_colors):
    kmer = draw(s.text('ACGT', min_size=kmer_size, max_size=kmer_size))
    coverage = tuple(
        draw(s.lists(s.integers(min_value=0), min_size=num_colors, max_size=num_colors)))
    edges = draw(s.lists(
        s.lists(s.booleans(), min_size=8, max_size=8),
        min_size=num_colors,
        max_size=num_colors))
    edges = tuple([tuple(edge_set) for edge_set in edges])
    return KmerRecord(tuple(kmer), coverage, edges)


@attr.s(slots=True)
class CortexGraphHeaderStub(object):
    kmer_size = attr.ib()
    kmer_container_size = attr.ib()
    num_colors = attr.ib()


class TestGraphBodyParser(object):
    @given(s.data(),
           s.integers(min_value=1),
           s.integers(min_value=0),
           s.integers(min_value=0, max_value=4))
    def test_parses_record(self, data, kmer_size, num_colors, n_kmers):
        # given
        kmer_container_size = max(1, ceil(kmer_size / 32))
        header = CortexGraphHeaderStub(kmer_size, kmer_container_size, num_colors)
        builder = CortexGraphBodyBuilder(kmer_container_size)

        expected_kmers = []
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        # when
        for kmer, expected_kmer in zip(
                kmer_generator_from_stream(BytesIO(builder.build()), header), expected_kmers):
            # then
            assert expected_kmer.kmer == kmer.kmer
            assert expected_kmer.coverage == kmer.coverage
            assert expected_kmer.edges == kmer.edges

    def test_parses_kmer_with_three_letters(self):
        kmer_container_size = 1
        kmer_size = 3
        num_colors = 0
        header = CortexGraphHeaderStub(kmer_size, kmer_container_size, num_colors)
        builder = CortexGraphBodyBuilder(kmer_container_size)

        expected_kmer = KmerRecord(tuple('AAC'), tuple(), tuple())
        builder.with_kmer_record(expected_kmer)

        kmer = next(kmer_generator_from_stream(BytesIO(builder.build()), header))

        assert expected_kmer.kmer == kmer.kmer
        assert expected_kmer.coverage == kmer.coverage
        assert expected_kmer.edges == kmer.edges
