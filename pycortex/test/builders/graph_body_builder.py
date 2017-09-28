import attr
from attr import Factory
from math import ceil
from bitstring import BitArray
from struct import pack

from hypothesis import strategies as s
from io import BytesIO

KMER_LETTER_TO_NUM = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
KMER_LETTER_TO_NUM_REVERSED_BITS = {'A': 0, 'C': 2, 'G': 1, 'T': 3}
KMER_LETTER_NUM_TO_BIT_REPR = ('00', '01', '10', '11')
KMER_CONTAINER_WORD_SIZE_IN_BITS = 64


def pack_edge_set(edge_set):
    edge_set = [str(flag) for flag in edge_set]
    edge_set = edge_set[:4] + list(reversed(edge_set[4:]))
    edge = BitArray('0b' + ''.join(edge_set))
    return edge.tobytes()


@attr.s(slots=True)
class KmerRecord(object):
    kmer = attr.ib()
    coverage = attr.ib()
    edges = attr.ib()
    kmer_container_size_in_uint64ts = attr.ib(1)

    def __attrs_post_init__(self):
        assert len(self.edges) == len(self.coverage)

    def to_bytestring(self):
        byte_string = self.pack_kmer()
        for color_coverage in self.coverage:
            byte_string += pack('I', color_coverage)
        for edge_set in self.edges:
            byte_string += pack_edge_set(edge_set)
        return byte_string

    def pack_kmer(self):
        kmer_string = ''.join(self.kmer).upper()
        bit_groups = []
        for letter_group_idx in range(ceil(len(kmer_string) / 4)):
            bit_group = []
            for letter_num in range(4):
                letter_idx = letter_num + letter_group_idx * 4
                if letter_idx >= len(kmer_string):
                    bit_group.insert(0, '00')
                else:
                    bit_group.append(KMER_LETTER_NUM_TO_BIT_REPR[
                                         KMER_LETTER_TO_NUM[kmer_string[letter_idx]]
                                     ])
            bit_groups.extend(bit_group)
        expected_packed_kmer_size_in_bit_groups = self.kmer_container_size_in_uint64ts * 32
        assert len(bit_groups) <= expected_packed_kmer_size_in_bit_groups
        missing_bit_groups = expected_packed_kmer_size_in_bit_groups - len(bit_groups)
        for _ in range(missing_bit_groups):
            bit_groups.append('00')
        bits = BitArray('0b' + ''.join(bit_groups))
        return bits.tobytes()


@attr.s(slots=True)
class CortexGraphBodyBuilder(object):
    sort_kmers = attr.ib(False)
    kmer_records = attr.ib(Factory(list), init=False)
    kmer_size = attr.ib(None)
    _kmer_container_size = attr.ib(None)

    @property
    def kmer_container_size(self):
        if self._kmer_container_size is None:
            if self.kmer_size is None:
                raise Exception("Need to call with_kmer_record to set kmer size first")
            self._kmer_container_size = max(1, ceil(self.kmer_size / 32))
        return self._kmer_container_size

    @kmer_container_size.setter
    def kmer_container_size(self, value):
        self._kmer_container_size = value

    def with_kmer_record(self, kmer):
        if self.kmer_size is None:
            self.kmer_size = len(kmer.kmer)
        assert len(kmer.kmer) == self.kmer_size

        coverage_per_color, edges_per_color = kmer.coverage, kmer.edges
        if isinstance(coverage_per_color, int):
            coverage_per_color = [coverage_per_color]
        if isinstance(edges_per_color, str):
            edges_per_color = [edges_per_color]
        assert self.kmer_container_size >= len(kmer.kmer) / 32
        self.kmer_records.append(
            KmerRecord(kmer.kmer, coverage_per_color, edges_per_color, self.kmer_container_size))
        return self

    def build(self):
        byte_strings = []
        if self.sort_kmers:
            self.kmer_records = sorted(self.kmer_records, key=lambda r: r.kmer)
        for record in self.kmer_records:
            byte_strings.append(record.to_bytestring())
        return BytesIO(b''.join(byte_strings))


@s.composite
def kmers(draw, kmer_size, num_colors):
    kmer = draw(s.text('ACGT', min_size=kmer_size, max_size=kmer_size))
    coverage = tuple(
        draw(s.lists(s.integers(min_value=0), min_size=num_colors, max_size=num_colors)))
    edges = draw(s.lists(
        s.lists(s.integers(min_value=0, max_value=1), min_size=8, max_size=8),
        min_size=num_colors,
        max_size=num_colors))
    edges = tuple([tuple(edge_set) for edge_set in edges])
    return KmerRecord(kmer, coverage, edges)
