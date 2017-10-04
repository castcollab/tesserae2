import attr
import struct
from bisect import bisect_left
from collections import Sequence
from io import SEEK_END
from struct import unpack
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
UINT8_T = 1
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
    record_size = attr.ib()

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

        kmer_container_size = header[2]
        if kmer_container_size <= 0:
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

        record_size = UINT64_T * kmer_container_size + (UINT32_T + UINT8_T) * num_colors
        return CortexGraphHeader(version=version, kmer_size=kmer_size,
                                 kmer_container_size=kmer_container_size,
                                 num_colors=num_colors,
                                 mean_read_lengths=mean_read_lengths,
                                 mean_total_sequence=mean_total_sequence,
                                 sample_names=sample_names,
                                 record_size=record_size)


def kmer_generator_from_stream(stream, cortex_header):
    record_size = cortex_header.kmer_container_size * UINT64_T + 5 * cortex_header.num_colors

    raw_record = stream.read(record_size)
    while raw_record != b'':
        yield CortexKmer(raw_record,
                         cortex_header.kmer_size,
                         cortex_header.num_colors,
                         cortex_header.kmer_container_size)
        raw_record = stream.read(record_size)


def edge_set_as_string(edge_set, is_revcomp=False):
    letters = []

    if is_revcomp:
        num_to_letter = list(reversed(NUM_TO_LETTER))
    else:
        num_to_letter = NUM_TO_LETTER

    for idx, edge in enumerate(edge_set):
        letter = num_to_letter[idx % 4]
        if idx < 4:
            letter = letter.lower()
        if edge:
            letters.append(letter)
        else:
            letters.append('.')

    if is_revcomp:
        incoming, outgoing = letters[:4], letters[4:]
        incoming, outgoing = list(reversed(incoming)), list(reversed(outgoing))
        letters = outgoing + incoming

    return ''.join(letters)


def cortex_kmer_as_cortex_jdk_print_string(kmer, alt_kmer_string=None):
    if kmer is None:
        revcomp_kmer = revcomp(alt_kmer_string)
        if revcomp_kmer > alt_kmer_string:
            revcomp_kmer = alt_kmer_string
        return '{}: {} missing'.format(revcomp_kmer, alt_kmer_string)
    if alt_kmer_string is not None and kmer.kmer != alt_kmer_string:
        is_revcomp = True
    else:
        is_revcomp = False

    edge_set_strings = [edge_set_as_string(edge_set, is_revcomp=is_revcomp) for edge_set in
                        kmer.edges]
    to_print = [str(kmer.kmer)]
    if alt_kmer_string is not None:
        to_print.append(': ' + alt_kmer_string)
    to_print.append(' ' + ' '.join(map(str, kmer.coverage)))
    to_print.append(' ' + ' '.join(edge_set_strings))
    return ''.join(to_print)


@attr.s(slots=True)
class CortexKmer(object):
    _raw_data = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    kmer_container_size_in_uint64ts = attr.ib(1)
    _kmer = attr.ib(None)
    _coverage = attr.ib(None)
    _edges = attr.ib(None)
    _kmer_vals_to_delete = attr.ib(init=False)

    def __attrs_post_init__(self):
        n_vals_left_over = self.kmer_size % 4
        n_vals_to_remove = 4 - n_vals_left_over
        if n_vals_to_remove > 0:
            self._kmer_vals_to_delete = (
                np.arange(0, n_vals_to_remove) + self.kmer_size - n_vals_left_over
            )

    @property
    def kmer(self):
        if self._kmer is None:
            kmer_as_uint64ts = np.frombuffer(
                self._raw_data[:self.kmer_container_size_in_uint64ts * 8],
                dtype='<u8')
            kmer_as_uint64ts_be = kmer_as_uint64ts.byteswap().newbyteorder()  # change to big endian
            kmer_as_properly_ordered_bits_right_aligned = np.unpackbits(
                np.frombuffer(kmer_as_uint64ts_be.tobytes(), dtype=np.uint8)
            )
            kmer = (
                kmer_as_properly_ordered_bits_right_aligned.reshape(-1, 2) * np.array([2, 1])
            ).sum(1)
            self._kmer = ''.join(NUM_TO_LETTER[num] for num in kmer[(len(kmer) - self.kmer_size):])
        return self._kmer

    @property
    def coverage(self):
        if self._coverage is None:
            start = self.kmer_container_size_in_uint64ts * UINT64_T
            coverage_raw = self._raw_data[start:(start + self.num_colors * UINT32_T)]
            fmt_string = ''.join(['I' for _ in range(self.num_colors)])
            self._coverage = unpack(fmt_string, coverage_raw)
        return self._coverage

    @property
    def edges(self):
        if self._edges is None:
            start = (
                self.kmer_container_size_in_uint64ts * UINT64_T + self.num_colors * UINT32_T
            )
            edge_bytes = list(self._raw_data[start:])
            edge_sets = np.unpackbits(np.array(edge_bytes, dtype=np.uint8)).reshape(-1, 4)
            edge_sets[1::2] = np.fliplr(edge_sets[1::2])
            self._edges = tuple(map(tuple, edge_sets.reshape(-1, 8)))
        return self._edges

    def get_edge_number_for_color(self, edge_num, color_num=0):
        return self.edges[color_num][edge_num]


class CortexKmerComparator(object):
    def __init__(self, *, kmer=None, kmer_object=None):
        self.kmer = kmer
        self.kmer_object = kmer_object
        if self.kmer is None:
            self.kmer = self.kmer_object.kmer

    def __eq__(self, other):
        return self.kmer == other.kmer

    def __lt__(self, other):
        return self.kmer < other.kmer

    def __repr__(self):
        return self.kmer


@attr.s(slots=True)
class CortexGraphStreamingParser(object):
    fh = attr.ib()
    header = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.header = CortexGraphHeader.from_stream(self.fh)

    def kmers(self):
        return kmer_generator_from_stream(self.fh, self.header)


@attr.s()
class CortexGraphKmerRecordSequence(Sequence):
    fh = attr.ib()
    cortex_header = attr.ib()
    body_start = attr.ib()
    n_records = attr.ib()
    record_size = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    kmer_size = attr.ib(init=False)
    kmer_container_size = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.record_size = self.cortex_header.record_size
        self.kmer_size = self.cortex_header.kmer_size
        self.num_colors = self.cortex_header.num_colors
        self.kmer_container_size = self.cortex_header.kmer_container_size

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError("Index must be of type int")
        if item >= self.n_records or item < 0:
            raise IndexError("Index ({}) is out of range".format(item))
        self.fh.seek(self.body_start + self.record_size * item)
        kmer_bytes = self.fh.read(self.record_size)
        return CortexKmerComparator(
            kmer_object=CortexKmer(
                kmer_bytes,
                kmer_size=self.kmer_size,
                num_colors=self.num_colors,
                kmer_container_size_in_uint64ts=self.kmer_container_size,
            )
        )

    def __len__(self):
        return self.n_records


def revcomp(kmer_string):
    return str(Seq(kmer_string, IUPAC.unambiguous_dna).reverse_complement())


class CortexGraphRandomAccessError(KeyError):
    """Raise this if a random access cortex graph parser could not find a kmer"""


@attr.s(slots=True)
class CortexGraphRandomAccessParser(object):
    fh = attr.ib()
    header = attr.ib(init=False)
    graph_sequence = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.fh.seekable()
        self.fh.seek(0)
        self.header = CortexGraphHeader.from_stream(self.fh)
        body_start_stream_position = self.fh.tell()

        self.fh.seek(0, SEEK_END)
        body_size = self.fh.tell() - body_start_stream_position
        if body_size % self.header.record_size != 0:
            raise RuntimeError(
                "Body size ({}) % Record size ({}) != 0".format(body_size, self.header.record_size))
        n_records = body_size // self.header.record_size
        self.graph_sequence = CortexGraphKmerRecordSequence(fh=self.fh,
                                                            body_start=body_start_stream_position,
                                                            cortex_header=self.header,
                                                            n_records=n_records)

    def get_kmer(self, kmer_string):
        kmer = CortexKmerComparator(kmer=kmer_string)
        try:
            kmer_comparator = index(self.graph_sequence, kmer, retrieve=True)
        except ValueError as e:
            raise CortexGraphRandomAccessError('Could not retrieve kmer: ' + kmer_string) from e

        return kmer_comparator.kmer_object

    def get_kmer_for_string(self, kmer_string):
        """Will compute the revcomp of kmer string before getting a kmer"""
        kmer_string_revcomp = revcomp(kmer_string)
        if kmer_string < kmer_string_revcomp:
            return self.get_kmer(kmer_string)
        else:
            return self.get_kmer(kmer_string_revcomp)


# copied from https://docs.python.org/3.6/library/bisect.html
def index(a, x, retrieve=False):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        val = a[i]
        if val == x:
            if retrieve:
                return val
            else:
                return i
    raise ValueError("Could not find '{}'".format(x))
