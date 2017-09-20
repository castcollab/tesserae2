import io
from struct import pack

import attr
from attr import Factory

"""
0	1	_Bool	Is tip-clipping applied?
1	1	_Bool	Have low-coverage unitigs been removed?
2	1	_Bool	Have low-coverage kmers been removed?
3	1	_Bool	Has this graph been cleaned against another graph?
4	4	uint32_t	Coverage threshold on unitigs
8	4	uint32_t	Coverage threshold on kmers
12	4	uint32_t	Length G of name of graph against which this graph has been cleaned
16	G	char[G]	Name of graph against which this graph has been cleaned
"""


@attr.s(slots=True)
class ColorInformationBlock(object):
    is_clipping_applied = attr.ib()
    are_low_coverage_unitigs_removed = attr.ib()
    are_low_coverage_kmers_removed = attr.ib()
    is_cleaned_against_other_graph = attr.ib()
    coverage_threshold_on_unitigs = attr.ib()
    coverage_threshold_on_kmers = attr.ib()
    length_of_cleaning_graph_name = attr.ib()
    cleaning_graph_name = attr.ib()

    def to_binary(self):
        string = pack('ccccIII', self.is_clipping_applied,
                      self.are_low_coverage_unitigs_removed,
                      self.are_low_coverage_kmers_removed,
                      self.is_cleaned_against_other_graph,
                      self.coverage_threshold_on_unitigs,
                      self.coverage_threshold_on_kmers,
                      self.length_of_cleaning_graph_name)
        return string + self.cleaning_graph_name


@attr.s(slots=True)
class CortexGraphBuilder(object):
    magic_word = attr.ib(b'CORTEX')
    version = attr.ib(6)
    kmer_size = attr.ib(1)
    kmer_bits = attr.ib(1)
    num_colors = attr.ib(1)
    error_rate = attr.ib(bytes(16))
    _mean_read_lengths = attr.ib(None)
    _total_sequence = attr.ib(None)
    _color_names = attr.ib(None)
    _color_information_blocks = attr.ib(Factory(list))

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
    def color_names(self):
        if self._color_names is None:
            self._color_names = ['sample_{}'.format(c).encode() for c in
                                 range(self.num_colors)]
        return self._color_names

    @property
    def color_information_blocks(self):
        if len(self._color_information_blocks) == 0:
            zero_bytes = [b'\x00' for _ in range(4)]
            self._color_information_blocks = [ColorInformationBlock(*zero_bytes, 0, 0, 0, b'')
                                              for _ in range(self.num_colors)]
        return self._color_information_blocks

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

    def with_color_names(self, color_names):
        self._color_names = color_names
        return self

    def with_color_information_block(self, color_information_block):
        self._color_information_blocks.append(color_information_block)
        return self

    def with_error_rate(self, error_rate):
        self.error_rate = error_rate
        return self

    def build(self):
        header_string = (
            self.magic_word +
            pack('IIII', self.version, self.kmer_size, self.kmer_bits, self.num_colors) +
            pack('{}I'.format(len(self.mean_read_lengths)), *self.mean_read_lengths) +
            pack('{}L'.format(len(self.total_sequence)), *self.total_sequence)
        )

        assert len(self.color_names) == len(self.color_information_blocks)
        assert len(self.color_names) == self.num_colors

        for color in self.color_names:
            header_string += pack('I', len(color))
            header_string += color

        header_string += self.error_rate

        for block in self.color_information_blocks:
            header_string += block.to_binary()

        header_string += self.magic_word

        return io.BytesIO(header_string)
