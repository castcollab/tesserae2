import struct

import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import strategies as s
from hypothesis.strategies import data, composite, binary, integers

from pycortex.graph_parser import CortexGraph, CortexGraphParserException
from pycortex.test.builders.cortex_graph_builder import CortexGraphBuilder, ColorInformationBlock

MAX_UINT = 2 ** (struct.calcsize('I') * 8) - 1
MAX_ULONG = 2 ** (struct.calcsize('L') * 8) - 1


@composite
def color_information_blocks(draw):
    bools = [draw(binary(min_size=1, max_size=1)) for _ in range(4)]
    uint32_ts = [draw(integers(min_value=0)) for _ in range(2)]
    name_size = draw(integers(min_value=0, max_value=3))
    name = draw(binary(min_size=name_size, max_size=name_size))

    return ColorInformationBlock(*bools, *uint32_ts, name_size, name)


class TestCortexGraphParsing(object):
    @given(s.binary())
    def test_raises_on_incorrect_magic_word(self, magic_word):
        assume(magic_word != b'CORTEX')

        fh = CortexGraphBuilder().with_magic_word(magic_word).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Saw magic word' in str(excinfo.value)

    @given(s.integers(min_value=0, max_value=MAX_UINT))
    def test_raises_on_incorrect_version(self, version):
        assume(version != 6)

        fh = CortexGraphBuilder().with_version(version).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Saw version' in str(excinfo.value)

    def test_raises_on_invalid_kmer_size(self):
        fh = CortexGraphBuilder().with_kmer_size(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Saw kmer size' in str(excinfo.value)

    def test_raises_on_invalid_kmer_bits(self):
        fh = CortexGraphBuilder().with_kmer_size(3).with_kmer_bits(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Saw kmer bits' in str(excinfo.value)

    def test_raises_on_invalid_num_colors(self):
        fh = CortexGraphBuilder().with_kmer_size(3).with_kmer_bits(1).with_num_colors(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Saw number of colors' in str(excinfo.value)

    @given(s.integers(min_value=1, max_value=10))
    def test_raises_when_concluding_magic_word_is_wrong(self, num_colors):
        fh = (CortexGraphBuilder()
              .with_num_colors(num_colors)
              .with_mean_read_lengths([0 for _ in range(num_colors + 1)])
              .build())

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraph.from_stream(fh)

        assert 'Concluding magic word' in str(excinfo.value)

    @given(data())
    def test_loads_entire_header_successfully(self, data):
        num_colors = data.draw(s.integers(min_value=1, max_value=3))
        kmer_size = data.draw(s.integers(min_value=1, max_value=100))
        kmer_bits = data.draw(s.integers(min_value=1, max_value=5))

        mean_read_lengths = data.draw(
            s.lists(elements=s.integers(min_value=0, max_value=MAX_UINT),
                    min_size=num_colors, max_size=num_colors))

        total_sequence = data.draw(
            s.lists(elements=s.integers(min_value=0, max_value=MAX_ULONG),
                    min_size=num_colors, max_size=num_colors))

        color_names = data.draw(s.lists(
            elements=s.binary(min_size=1, max_size=3), min_size=num_colors, max_size=num_colors)
        )

        cgb = (CortexGraphBuilder()
               .with_kmer_size(kmer_size)
               .with_kmer_bits(kmer_bits)
               .with_num_colors(num_colors)
               .with_mean_read_lengths(mean_read_lengths)
               .with_total_sequence(total_sequence)
               .with_color_names(color_names)
               .with_error_rate(data.draw(s.binary(min_size=16, max_size=16))))

        for color_name in color_names:
            cgb.with_color_information_block(data.draw(color_information_blocks()))

        fh = cgb.build()

        cp = CortexGraph.from_stream(fh)
