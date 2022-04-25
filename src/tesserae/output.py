from __future__ import annotations
import sys
from typing import Optional

import pysam

from tesserae.utils import RefInterval, calc_num_mismatches


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
        segment.set_tag('NM', nm)
        segment.set_tag('LL', float(log_likelihood), 'f')

        self.bam.write(segment)
