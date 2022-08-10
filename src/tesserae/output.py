from __future__ import annotations
import functools
import sys
from typing import Optional, TYPE_CHECKING

import pysam

from tesserae.utils import RefInterval, calc_num_mismatches

if TYPE_CHECKING:
    from skbio import Sequence


def pprint_alignment(bam: pysam.AlignmentFile, refs: list[Sequence], file=None):
    print_file = functools.partial(print, file=file)

    ref_seq_by_id = {r.metadata['id']: str(r).upper() for r in refs}
    ivals = list(sorted(bam.fetch(until_eof=True), key=lambda aln: aln.get_tag('IV') if aln.has_tag('IV') else 0))

    # Prepare lists
    qry_alignment = []
    aln_symbols = []
    ref_alignments = []

    # Parse all aligned intervals
    curr_aln_pos = 0
    query_id = ""
    query_len = 0
    query_log_likelihood = float('-inf')
    for aln in ivals:
        query_id = aln.query_name
        query_len += aln.query_length
        query_log_likelihood = aln.get_tag('LL')

        ref_alignment = []
        curr_ref_pos = aln.reference_start
        curr_ref_id = aln.reference_name
        for qry_pos, ref_pos in aln.get_aligned_pairs():
            if qry_pos is not None and ref_pos is not None:
                qry_alignment.append(aln.query_sequence[qry_pos])
                ref_alignment.append(ref_seq_by_id[curr_ref_id][ref_pos])
                aln_symbols.append("|" if qry_alignment[-1] == ref_alignment[-1] else " ")
            elif qry_pos is not None and ref_pos is None:
                qry_alignment.append(aln.query_sequence[qry_pos])
                ref_alignment.append("-")
                aln_symbols.append("^")
            elif qry_pos is None and ref_pos is not None:
                qry_alignment.append("-")
                ref_alignment.append(ref_seq_by_id[curr_ref_id][ref_pos])
                aln_symbols.append("~")

        # Extend the reference sequence a bit, but mask it with lowercase to indicate those are unaligned
        ref_seq_extra = ref_seq_by_id[curr_ref_id][curr_ref_pos:curr_ref_pos + 25].lower()
        ref_alignment_str = "".join(ref_alignment)
        ref_alignments.append((curr_aln_pos, ref_alignment_str + ref_seq_extra))
        curr_aln_pos += len(ref_alignment_str)

    # Find longest sequence ID
    max_id_length = max(len(k) for k in ref_seq_by_id.keys())
    max_id_length = max(max_id_length, len(query_id))

    # Find number of characters required to display the interval positions
    start_pos_len = max(len(str(aln.reference_start)) for aln in ivals) if ivals else 1
    end_pos_len = max(len(str(aln.reference_end)) for aln in ivals) if ivals else 1
    end_pos_len = max(end_pos_len, len(str(query_len)))
    pos_len = max(start_pos_len, end_pos_len)

    qry_start_pos = 0
    qry_header = (f"{query_id:>{max_id_length}} [{qry_start_pos:>{pos_len}}"
                  f"-{query_len:<{pos_len}})    ")
    print_file(qry_header, "".join(qry_alignment))
    print_file(" " * len(qry_header), "".join(aln_symbols))

    for aln, (aln_start, ref_aln) in zip(ivals, ref_alignments):
        header = (f"{aln.reference_name:>{max_id_length}} [{aln.reference_start:>{pos_len}}"
                  f"-{aln.reference_end:<{pos_len}})    ")
        print_file(header, (" " * aln_start) + ref_aln)

    print_file("Log-likelihood:", query_log_likelihood)
    print_file()


class SamWriter:
    def __init__(self, ref_contigs, fname, mode, program):
        self.ref_contigs = ref_contigs
        self.bam_header = {
            'HD': {'VN': '1.6'},
            'SQ': [],
            'PG': [program]
        }

        self.ref_name_to_id = {}

        for i, r in enumerate(ref_contigs):
            self.bam_header['SQ'].append({
                'SN': r.metadata['id'],
                'LN': len(r)
            })

            self.ref_name_to_id[r.metadata['id']] = i

        self.fname = fname
        self.mode = mode

        self.bam: Optional[pysam.AlignmentFile] = None

    def __enter__(self):
        self.bam = pysam.AlignmentFile(self.fname, self.mode, header=self.bam_header)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.flush()
        self.bam.close()
        self.bam = None

        return self

    def write_alignment(self, query_id: str, query: str, interval: RefInterval, ival_num: int, log_likelihood: float):
        if not self.bam:
            raise IOError("BAM file not yet open! Use this class using a 'with' statement.")

        segment = pysam.AlignedSegment()

        segment.query_name = query_id
        segment.reference_id = self.ref_name_to_id[interval.ref_id]
        segment.reference_start = interval.ref_start if interval.ref_start else 0
        segment.cigarstring = interval.cigar_str()

        if ival_num > 0:
            segment.flag |= pysam.FSUPPLEMENTARY

        query_sub_seq = query[interval.qry_start:interval.qry_end]
        segment.query_sequence = query_sub_seq
        segment.mapping_quality = 60

        ref_seq = str(self.ref_contigs[segment.reference_id])[interval.ref_start:interval.ref_end]
        nm = calc_num_mismatches(query_sub_seq, ref_seq, interval.cigar_ops)
        snps = calc_num_mismatches(query_sub_seq, ref_seq, interval.cigar_ops, ignore_indels=True)
        segment.set_tag('NM', nm)
        segment.set_tag('NS', snps)
        segment.set_tag('LL', float(log_likelihood), 'f')
        segment.set_tag('IV', ival_num)

        self.bam.write(segment)
