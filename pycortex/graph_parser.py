import struct
from struct import unpack

CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6


class CortexGraphParserException(Exception):
    pass


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
