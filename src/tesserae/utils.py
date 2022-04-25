import re
import bz2
import gzip
import lzma
from pathlib import Path
from dataclasses import dataclass, field
from typing import Union, TextIO, BinaryIO, Iterable


def open_compressed(filename, *args, **kwargs) -> Union[TextIO, BinaryIO]:
    if not isinstance(filename, Path):
        filename = Path(filename)

    fopen = {
        '.gz': gzip.open,
        '.bz2': bz2.open,
        '.lz': lzma.open
    }.get(filename.suffix, open)

    return fopen(filename, *args, **kwargs)


@dataclass
class RefInterval:
    ref_id: str = None
    ref_start: int = None
    ref_end: int = None
    qry_start: int = None
    qry_end: int = None
    cigar_ops: list[tuple[int, str]] = field(default_factory=list[tuple[int, str]])

    def cigar_str(self) -> str:
        return "".join(f"{c}{op}" for c, op in self.cigar_ops)

    def alignment_length(self) -> int:
        return sum(c for c, op in self.cigar_ops)

    def __str__(self):
        return (f"RefInterval(ref_id={self.ref_id}, ref_start={self.ref_start}, ref_end={self.ref_end}, "
                f"qry_start={self.qry_start}, qry_end={self.qry_end}, cigar={self.cigar_str()})")

    def __repr__(self):
        return self.__str__()


cigar_op_re = re.compile(r"(\d+[MID])")


def cigar_str_to_tuples(cigar_str: str) -> Iterable[tuple[int, str]]:
    for match in cigar_op_re.finditer(cigar_str):
        op = match.group(1)
        yield int(op[:-1]), op[-1]


def calc_num_mismatches(query_seq: str, ref_seq: str, cigar_ops: list[tuple[int, str]], ignore_indels=False) -> int:
    nm = 0

    query_pos = 0
    ref_contig_pos = 0
    for i, (count, op) in enumerate(cigar_ops):
        if op == "M":
            nm += sum(query_seq[query_pos + j] != ref_seq[ref_contig_pos + j] for j in range(count))

            query_pos += count
            ref_contig_pos += count
        elif op == "I":
            if 0 < i < len(cigar_ops)-1 and not ignore_indels:
                nm += count

            query_pos += count
        elif op == "D" and not ignore_indels:
            nm += count
            ref_contig_pos += count
        elif op not in {"M", "I", "D"}:
            raise ValueError(f"Invalid CIGAR op {op}")

    return nm
