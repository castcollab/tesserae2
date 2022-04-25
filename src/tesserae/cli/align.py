from __future__ import annotations

import logging
import argparse
import sys
from pathlib import Path

import skbio  # type: ignore

from tesserae import TesseraeFactory, open_compressed, SamWriter, __version__
from tesserae.cli.registry import Subcommand

logger = logging.getLogger(__name__)


OUTPUT_TYPES = {
    "sam": "w",
    "bam": "wb"
}


class AlignSubcommand(Subcommand):
    """
    Align query sequences to a panel of reference sequences.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'ref_panel', type=Path,
            help="Path to a FASTA file with the reference sequences. Can be compressed with gzip or bzip2."
        )

        subparser.add_argument(
            'query_seqs', type=Path,
            help="Path to a FASTA file with the query sequences to align. Can be compressed with gzip or bzip2."
        )

        subparser.add_argument(
            '-O', '--output-type', default="sam", choices=OUTPUT_TYPES.keys(),
            help="Output type, options include sam or bam output."
        )

        subparser.add_argument(
            '-o', '--output', default="-",
            help="Output filename. Defaults to stdout."
        )

    def __call__(self, ref_panel: Path, query_seqs: Path, output_type: str, output: str,  # type: ignore
                 *args, **kwargs):  # type: ignore
        with open_compressed(ref_panel) as f:
            ref_contigs = list(skbio.io.read(f, "fasta"))

        # TODO: configure HMM params
        logger.info("Configuring HMM...")
        factory = TesseraeFactory()
        hmm = factory.build(ref_contigs)

        output_mode = OUTPUT_TYPES[output_type]

        program = {
            'ID': 'tesserae0',
            'PN': 'tesserae',
            'VN': __version__,
            'CL': " ".join(sys.argv)
        }

        with open_compressed(query_seqs) as f, SamWriter(ref_contigs, output, output_mode, program) as w:
            logger.info("Aligning query contigs...")
            for query in skbio.io.read(f, "fasta"):
                query_str = str(query).upper()
                result = hmm.align(query_str.encode('ascii'))

                for i, interval in enumerate(result.get_aligned_ref_intervals()):
                    w.write_alignment(query.metadata['id'], query_str, interval, i, result.get_log_likelihood())

        logger.info("Done.")
