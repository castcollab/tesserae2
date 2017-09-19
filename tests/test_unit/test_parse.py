import ctypes
from struct import pack, unpack

import io
import pytest
import struct

from hypothesis import assume
from hypothesis import given
from hypothesis import strategies as s
from hypothesis.strategies import data


class CortexGraphParserException(Exception):
    pass


CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
MAX_UINT = 2 ** (struct.calcsize('I') * 8) - 1
MAX_ULONG = 2 ** (struct.calcsize('L') * 8) - 1


class CortexGraphParser(object):
    def __init__(self, fh):
        self.fh = fh

        self._read_header()

    def _read_header(self):
        magic_word = unpack('cccccc', self.fh.read(6))

        if magic_word != CORTEX_MAGIC_WORD:
            raise CortexGraphParserException(
                "Saw magic word {} but was expecting {}".format(magic_word, CORTEX_MAGIC_WORD))

        header = unpack('IIII', self.fh.read(16))

        self.version = header[0]
        if self.version != CORTEX_VERSION:
            raise CortexGraphParserException(
                "Saw version {} but was expecting {}".format(self.version, CORTEX_VERSION)
            )

        self.kmer_size = header[1]
        if self.kmer_size <= 0:
            raise CortexGraphParserException(
                "Saw kmer size {} but was expecting value > 0".format(self.kmer_size)
            )

        self.kmer_bits = header[2]
        if self.kmer_bits <= 0:
            raise CortexGraphParserException(
                "Saw kmer bits {} but was expecting value > 0".format(self.kmer_size)
            )

        self.num_colors = header[3]
        if self.num_colors <= 0:
            raise CortexGraphParserException(
                "Saw number of colors {} but was expecting value > 0".format(self.kmer_size)
            )

        self.mean_read_lengths = unpack(
            '{}I'.format(self.num_colors), self.fh.read(struct.calcsize('I') * self.num_colors)
        )

        self.mean_total_sequence = unpack(
            '{}L'.format(self.num_colors), self.fh.read(struct.calcsize('L') * self.num_colors)
        )

        self.sample_names = []
        for _ in range(self.num_colors):
            rv = self.fh.read(struct.calcsize('I'))
            snlength = unpack('I', rv)[0]
            sample_name = unpack('{}c'.format(snlength), self.fh.read(snlength))
            self.sample_names.append(sample_name)

        concluding_magic_word = unpack('cccccc', self.fh.read(6))

        if concluding_magic_word != magic_word:
            raise CortexGraphParserException(
                "Concluding magic word {} != starting magic word {}".format(concluding_magic_word,
                                                                            magic_word))


class CortexGraphBuilder(object):
    def __init__(self):
        self.magic_word = b'CORTEX'
        self.version = 6
        self.kmer_size = 1
        self.kmer_bits = 1
        self.num_colors = 1
        self._mean_read_lengths = None
        self._total_sequence = None
        self._sample_names = None

    @property
    def mean_read_lengths(self):
        if self._mean_read_lengths is None:
            self._mean_read_lengths = [0 for _ in range(self.num_colors)]

        return self._mean_read_lengths

    @property
    def total_sequence(self):
        if self._total_sequence is None:
            self._total_sequence = [0 for _ in range(self.num_colors)]

        return self._total_sequence

    @property
    def sample_names(self):
        if self._sample_names is None:
            self._sample_names = ['sample_{}'.format(c).encode() for c in range(self.num_colors)]

        return self._sample_names

    def with_magic_word(self, magic_word):
        self.magic_word = magic_word
        return self

    def with_version(self, version):
        self.version = version
        return self

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_kmer_bits(self, kmer_bits):
        self.kmer_bits = kmer_bits
        return self

    def with_num_colors(self, num_colors):
        self.num_colors = num_colors
        return self

    def with_mean_read_lengths(self, mean_read_lengths):
        self._mean_read_lengths = mean_read_lengths
        return self

    def with_total_sequence(self, total_sequence):
        self._total_sequence = total_sequence
        return self

    def with_sample_names(self, sample_names):
        self._sample_names = sample_names
        return self

    def build(self):
        header_string = (
            self.magic_word +
            pack('IIII', self.version, self.kmer_size, self.kmer_bits, self.num_colors) +
            pack('{}I'.format(len(self.mean_read_lengths)), *self.mean_read_lengths) +
            pack('{}L'.format(len(self.total_sequence)), *self.total_sequence)
        )

        for sample in self.sample_names:
            header_string += pack('I', len(sample))
            header_string += sample

        header_string += self.magic_word

        return io.BytesIO(header_string)


class TestCortexGraphParsing(object):
    @given(s.binary())
    def test_raises_on_incorrect_magic_word(self, magic_word):
        assume(magic_word != b'CORTEX')

        fh = CortexGraphBuilder().with_magic_word(magic_word).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

        assert 'Saw magic word' in str(excinfo.value)

    @given(s.integers(min_value=0, max_value=MAX_UINT))
    def test_raises_on_incorrect_version(self, version):
        assume(version != 6)

        fh = CortexGraphBuilder().with_version(version).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

        assert 'Saw version' in str(excinfo.value)

    def test_raises_on_invalid_kmer_size(self):
        fh = CortexGraphBuilder().with_kmer_size(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

        assert 'Saw kmer size' in str(excinfo.value)

    def test_raises_on_invalid_kmer_bits(self):
        fh = CortexGraphBuilder().with_kmer_size(3).with_kmer_bits(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

        assert 'Saw kmer bits' in str(excinfo.value)

    def test_raises_on_invalid_num_colors(self):
        fh = CortexGraphBuilder().with_kmer_size(3).with_kmer_bits(1).with_num_colors(0).build()

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

        assert 'Saw number of colors' in str(excinfo.value)

    @given(s.integers(min_value=1, max_value=10))
    def test_raises_when_concluding_magic_word_is_wrong(self, num_colors):
        fh = (CortexGraphBuilder()
              .with_num_colors(num_colors)
              .with_mean_read_lengths([0 for _ in range(num_colors + 1)])
              .build())

        with pytest.raises(CortexGraphParserException) as excinfo:
            cp = CortexGraphParser(fh)

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

        sample_names = data.draw(s.lists(
            elements=s.binary(min_size=1, max_size=256), min_size=num_colors, max_size=num_colors)
        )

        fh = (CortexGraphBuilder()
              .with_kmer_size(kmer_size)
              .with_kmer_bits(kmer_bits)
              .with_num_colors(num_colors)
              .with_mean_read_lengths(mean_read_lengths)
              .with_total_sequence(total_sequence)
              .with_sample_names(sample_names)
              .build())

        cp = CortexGraphParser(fh)
