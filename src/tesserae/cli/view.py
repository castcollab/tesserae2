from __future__ import annotations

import os
import contextlib
import logging
import argparse
import sys
from pathlib import Path

import pysam
import skbio  # type: ignore

from tesserae import open_compressed, pprint_alignment
from tesserae.cli.registry import Subcommand

logger = logging.getLogger(__name__)


class ViewSubcommand(Subcommand):
    """
    Pretty print an alignment in text form.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'ref_panel', type=Path,
            help="Path to a FASTA with the reference contigs used for alignment."
        )

        subparser.add_argument(
            'alignment', default="-",
            help="Path to SAM/BAM file with the alignment. Defaults to stdin."
        )

        subparser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help="Output file. Defaults to stdout."
        )

    def __call__(self, ref_panel, alignment, output, *args, **kwargs):  # type: ignore
        with open_compressed(ref_panel) as f:
            ref_contigs = list(skbio.io.read(f, "fasta"))

        with open(os.devnull) as null:
            with contextlib.redirect_stderr(null):
                bam = pysam.AlignmentFile(alignment)

        pprint_alignment(bam, ref_contigs, output)
