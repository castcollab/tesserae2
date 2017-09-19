import io
from struct import pack


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
