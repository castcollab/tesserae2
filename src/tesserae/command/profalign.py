import argparse
import contextlib
import logging
import math
import os
import pathlib
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import pysam

from tesserae.nucleotide_sequence import NucleotideSequence

from .. import profmodel

LOGGER = logging.getLogger(__name__)


def profalign(args) -> None:
    """Main CLI call for the 'profalign' sub-command."""

    aligner = profmodel.Tesserae2()
    LOGGER.info("Aligning with Tesserae2...")
    start_time = time.time()
    aligner.align_from_fastx(args.query, args.sources)
    end_time = time.time()
    LOGGER.info("Alignment took %fs", end_time - start_time)
    # dump_results_to_log(target_alignment_results)


def main(raw_args):

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        prog="tesserae align",
        description="Reads two input fasta files - a reads file and a panel "
        "of source sequences.  The query is aligned to all source sequences "
        "simultaneously, allowing for recombination and trimming of irrelevant "
        "flanking sequences.  Output in SAM file format on stdout (and "
        "optionally to a file).",
        usage="align query and source sequences",
    )

    profalign_required_args = parser.add_argument_group("Required Arguments")
    profalign_required_args.add_argument(
        "-q", "--query", help="Query FASTA", required=True
    )
    profalign_required_args.add_argument(
        "-s", "--sources", help="Source panel FASTA", required=True
    )

    parser.add_argument(
        "-o",
        "--bamout",
        help="Output BAM file in which to store alignment results.",
        required=False,
    )

    # Parse args
    args = parser.parse_args(args=raw_args)

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(sys.argv))
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our sub-command:
    profalign(args)

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)
