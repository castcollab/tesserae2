import struct
from struct import unpack

import attr

CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
UINT32_T = 4
UINT64_T = 8
BITS_IN_BYTE = 8
NUM_TO_LETTER = ['A', 'C', 'G', 'T']


class CortexGraphParserException(Exception):
    pass


@attr.s(slots=True)
class CortexGraphHeader(object):
    version = attr.ib()
    kmer_size = attr.ib()
    kmer_container_size = attr.ib()
    num_colors = attr.ib()
    mean_read_lengths = attr.ib()
    mean_total_sequence = attr.ib()
    sample_names = attr.ib()

    @classmethod
    def from_stream(cls, fh):
        magic_word = unpack('cccccc', fh.read(6))

        if magic_word != CORTEX_MAGIC_WORD:
            raise CortexGraphParserException(
                "Saw magic word {} but was expecting {}".format(magic_word, CORTEX_MAGIC_WORD))

        header = unpack('IIII', fh.read(16))

        version = header[0]
        if version != CORTEX_VERSION:
            raise CortexGraphParserException(
                "Saw version {} but was expecting {}".format(version, CORTEX_VERSION)
            )

        kmer_size = header[1]
        if kmer_size <= 0:
            raise CortexGraphParserException(
                "Saw kmer size {} but was expecting value > 0".format(kmer_size)
            )

        kmer_bytes = header[2]
        if kmer_bytes <= 0:
            raise CortexGraphParserException(
                "Saw kmer bits {} but was expecting value > 0".format(kmer_size)
            )

        num_colors = header[3]
        if num_colors <= 0:
            raise CortexGraphParserException(
                "Saw number of colors {} but was expecting value > 0".format(kmer_size)
            )

        mean_read_lengths = unpack(
            '{}I'.format(num_colors), fh.read(struct.calcsize('I') * num_colors)
        )

        mean_total_sequence = unpack(
            '{}L'.format(num_colors), fh.read(struct.calcsize('L') * num_colors)
        )

        sample_names = []
        for _ in range(num_colors):
            sample_name_length_string = fh.read(struct.calcsize('I'))
            snlength = unpack('I', sample_name_length_string)[0]
            sample_name = unpack('{}c'.format(snlength), fh.read(snlength))
            sample_names.append(b''.join(sample_name))
        sample_names = tuple(sample_names)

        error_rate = unpack('16c', fh.read(16))

        for _ in range(num_colors):
            color_info_block_string = fh.read(4 + 3 * struct.calcsize('I'))
            color_info_block = unpack('ccccIII', color_info_block_string)
            fh.read(color_info_block[6])

        concluding_magic_word = unpack('cccccc', fh.read(6))

        if concluding_magic_word != magic_word:
            raise CortexGraphParserException(
                "Concluding magic word {} != starting magic word {}".format(concluding_magic_word,
                                                                            magic_word))
        return CortexGraphHeader(version=version, kmer_size=kmer_size,
                                 kmer_container_size=kmer_bytes,
                                 num_colors=num_colors,
                                 mean_read_lengths=mean_read_lengths,
                                 mean_total_sequence=mean_total_sequence,
                                 sample_names=sample_names)


def kmer_generator_from_stream(stream, cortex_header):
    record_size = cortex_header.kmer_container_size * UINT64_T + 5 * cortex_header.num_colors

    while True:
        raw_record = stream.read(record_size)
        if raw_record == '':
            return
        yield CortexKmer(raw_record,
                         cortex_header.kmer_size,
                         cortex_header.num_colors,
                         cortex_header.kmer_container_size)


@attr.s(slots=True)
class CortexKmer(object):
    _raw_data = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    kmer_container_size_in_uint64ts = attr.ib(1)
    _kmer = attr.ib(None)
    _coverage = attr.ib(None)
    _edges = attr.ib(None)

    @property
    def kmer(self):
        if self._kmer is None:
            letters = []
            for index, byte in enumerate(self._raw_data[:self.kmer_container_size_in_uint64ts * 8]):
                if index * 4 > self.kmer_size:
                    break
                letters.append(NUM_TO_LETTER[(byte & 0b11000000) >> 6])
                letters.append(NUM_TO_LETTER[(byte & 0b110000) >> 4])
                letters.append(NUM_TO_LETTER[(byte & 0b1100) >> 2])
                letters.append(NUM_TO_LETTER[byte & 0b11])
            self._kmer = tuple(letters[:self.kmer_size])
        return self._kmer

    @property
    def coverage(self):
        if self._coverage is None:
            start = self.kmer_container_size_in_uint64ts * UINT64_T
            coverage_raw = self._raw_data[start:(start + self.num_colors * UINT32_T)]
            fmt_string = ''.join(['I' for _ in range(self.num_colors)])
            self._coverage = unpack(fmt_string, coverage_raw)
        return self._coverage

    def get_coverage_for_color(self, color_num):
        assert self.num_colors > color_num
        return self.coverage[color_num]

    @property
    def edges(self):
        if self._edges is None:
            start = (
                self.kmer_container_size_in_uint64ts * UINT64_T + self.num_colors * UINT32_T
            )
            edge_bytes = self._raw_data[start:]
            edges_by_color = []
            for edge_color_num in range(self.num_colors):
                edges = []
                for bit_idx in range(BITS_IN_BYTE):
                    edges.append((edge_bytes[edge_color_num] >> bit_idx) & 1 == 1)
                edges_by_color.append(tuple(edges[::-1]))
            self._edges = tuple(edges_by_color)
        return self._edges

    def get_edge_number_for_color(self, edge_num, color_num=0):
        return self.edges[color_num][edge_num]


@attr.s(slots=True)
class CortexGraph(object):
    header = attr.ib()

    @classmethod
    def from_stream(self, fh):
        header = CortexGraphHeader.from_stream(fh)
        return CortexGraph(header=header)
