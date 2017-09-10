import ctypes
from struct import pack, unpack

import io
import pytest
import struct

from hypothesis import assume
from hypothesis import example
from hypothesis import given
from hypothesis import strategies as s


class CortexGraphParserException(Exception):
    pass


CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
MAX_UINT = 2 ** (struct.calcsize('I') * 8) - 1


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
                "Saw kmer size {} but was expecting value >= 0".format(self.kmer_size)
            )

        self.kmer_bits = header[2]
        if self.kmer_bits <= 0:
            raise CortexGraphParserException(
                "Saw kmer bits {} but was expecting value >= 0".format(self.kmer_size)
            )

        self.num_colors = header[3]
        if self.num_colors <= 0:
            raise CortexGraphParserException(
                "Saw number of colors {} but was expecting value >= 0".format(self.kmer_size)
            )


class CortexGraphBuilder(object):
    def __init__(self):
        self.magic_word = b'CORTEX'
        self.version = 6
        self.kmer_size = 0
        self.kmer_bits = 0
        self.num_colors = 0

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

    def build(self):
        return io.BytesIO(
            self.magic_word +
            pack('IIII', self.version, self.kmer_size, self.kmer_bits, self.num_colors)
        )


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

    @given(kmer_size=s.integers(min_value=1, max_value=MAX_UINT),
           kmer_bits=s.integers(min_value=1, max_value=MAX_UINT),
           num_colors=s.integers(min_value=1, max_value=MAX_UINT))
    def test_loads_header_successfully(self, kmer_size, kmer_bits, num_colors):
        fh = (
            CortexGraphBuilder()
            .with_kmer_size(kmer_size)
            .with_kmer_bits(kmer_bits)
            .with_num_colors(num_colors)
            .build()
        )

        cp = CortexGraphParser(fh)
