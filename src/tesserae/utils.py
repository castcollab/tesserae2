from __future__ import annotations

import re
import bz2
import gzip
import lzma
from pathlib import Path
from dataclasses import dataclass, field
from typing import Union, TextIO, BinaryIO, Iterable, Optional


def open_compressed(filename, *args, **kwargs) -> Union[TextIO, BinaryIO]:
    if not isinstance(filename, Path):
        filename = Path(filename)

    fopen = {  # type: ignore
        '.gz': gzip.open,
        '.bz2': bz2.open,
        '.lz': lzma.open
    }.get(filename.suffix, open)

    return fopen(filename, *args, **kwargs)  # type: ignore


def parse_num_suffix(num, bytes=False):
    """
    Parse a string containing a number with a possible suffix like 'M' or 'G',
    and multiply accordingly.

    >>> parse_num_suffix('10M')
    10000000
    >>> parse_num_suffix('5G')
    5000000000
    >>> parse_num_suffix('3k')
    3000
    >>> parse_num_suffix('500')
    500

    Parameters
    ----------
    num : str
        The number with possible suffix to parse

    Returns
    -------
    int
        An integer multiplied accordingly to the suffix
    """

    if not num:
        return None

    base = 1024 if bytes else 1000

    suffixes = ['K', 'M', 'G', 'T', 'P']
    expo = list(range(1, len(suffixes)+1))

    suffix_expo = dict(zip(suffixes, expo))

    if not num[-1].isalpha():
        return int(num)

    suffix = num[-1].upper()
    if suffix not in suffixes:
        raise ValueError(
            "'{}' is not a valid number. Supported suffixes: {}".format(
                num, ", ".join(iter(suffixes.keys()))
            ))

    return int(num[:-1]) * (base ** suffix_expo[suffix])


@dataclass
class RefInterval:
    ref_id: Optional[str] = None
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None
    qry_start: Optional[int] = None
    qry_end: Optional[int] = None
    cigar_ops: list[tuple[int, str]] = field(default_factory=list)

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
            if 0 < i < len(cigar_ops) - 1 and not ignore_indels:
                nm += count

            query_pos += count
        elif op == "D" and not ignore_indels:
            nm += count
            ref_contig_pos += count
        elif op not in {"M", "I", "D"}:
            raise ValueError(f"Invalid CIGAR op {op}")

    return nm
