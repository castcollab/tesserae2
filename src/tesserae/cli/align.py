from __future__ import annotations

import logging
import argparse
import sys
from pathlib import Path

if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib

import numpy
import skbio  # type: ignore

from tesserae import TesseraeFactory, open_compressed, parse_num_suffix, SamWriter, __version__
from tesserae.cli.registry import Subcommand

logger = logging.getLogger(__name__)


OUTPUT_TYPES = {
    "sam": "w",
    "bam": "wb"
}


def add_hmm_config_options(subparser: argparse.ArgumentParser):
    subparser.add_argument(
        '-c', '--hmm-config', type=argparse.FileType('rb'), default=None, required=False,
        help="Configure HMM transmission and emission probabilities using the specified TOML config file."
    )

    hmm_group = subparser.add_argument_group("Hidden Markov Model configuration",  # noqa: F841
                                             "Any configuration specified on the command line overwrites "
                                             "configuration specified in a config file.")
    pi_group = subparser.add_argument_group("π_m and π_i, the probabilities of starting in a match or insert "
                                            "state respectively. Should sum to 1, so by specifying π_m, "
                                            "π_i is implicitly inferred")
    pi_group.add_argument(
        '--PM', type=float, required=False, metavar="π_m",
        help="Probability of starting in a match state (π_m)"
    )

    transition_group = subparser.add_argument_group("Gap open, gap extension, reference jump and termination "
                                                    "probabilities")
    transition_group.add_argument(
        '-D', '--gap-open', type=float, required=False, metavar="δ",
        help="Gap open probability (δ)."
    )
    transition_group.add_argument(
        '-E', '--gap-extension', type=float, required=False, metavar="ε",
        help="Gap extension probability (ε)"
    )
    transition_group.add_argument(
        '-R', '--reference-jump', type=float, required=False, metavar="ρ",
        help="Reference jump probability (ρ)"
    )
    transition_group.add_argument(
        '-T', '--termination', type=float, required=False, metavar="τ",
        help="Termination probability (τ)"
    )

    emission_group = subparser.add_argument_group(
        "Emission probabilities",
        "Match + TS + 2*TV probabilities should sum to one. By setting the match and TS probabilities, "
        "TV probability is implicitly inferred. If configuring emission probabilities, you'll have to specify both "
        "match and TS probabilities, you can't specify just one of them"
    )

    emission_group.add_argument(
        '-M', '--emission-match', type=float, required=False,
        help="Probability of emitting a matching pair of nucleotides."
    )
    emission_group.add_argument(
        '-S', '--emission-ts', type=float, required=False,
        help="Probability of emitting a transition, i.e., purine-purine or pyrimidine-pyrimidine conversion "
             "(A<=>G or C<=>T)."
    )
    emission_group.add_argument(
        '-I', '--emission-indel', type=float, nargs=4, required=False, metavar=('A', 'C', 'G', 'T'),
        help="Nucleotide distribution in an insert/deletion state, should specify four probabilities (that sum to "
             "one), representing the emission probability of A, C, G and T in an indel state."
    )


def configure_hmm_factory(hmm_params):
    factory = TesseraeFactory()

    if hmm_params.get('PM') is not None:
        factory.pi_m = hmm_params['PM']

    if hmm_params.get('gap_open') is not None:
        factory.delta = hmm_params['gap_open']

    if hmm_params.get('gap_extension') is not None:
        factory.eps = hmm_params['gap_extension']

    if hmm_params.get('reference_jump') is not None:
        factory.rho = hmm_params['reference_jump']

    if hmm_params.get('termination') is not None:
        factory.tau = hmm_params['termination']

    emission_match = hmm_params.get('emission_match')
    emission_ts = hmm_params.get('emission_ts')
    if emission_match is not None and emission_ts is None:
        logger.error("When configuring match emission probabilities, both match and transition probabilities are "
                     "required (options emission_match and emission_ts).")
        return 1
    elif emission_match is None and emission_ts is not None:
        logger.error("When configuring match emission probabilities, both match and transition probabilities are "
                     "required (options emission_match and emission_ts).")
        return 1
    elif emission_match is not None and emission_ts is not None:
        factory.configure_emissions(emission_match, emission_ts)

    if hmm_params.get('emission_indel') is not None:
        e_indel = numpy.array(hmm_params['emission_indel'], dtype=numpy.float32)
        if e_indel.sum() != 1.0:
            logger.error("Indel emission distribution doesn't sum to one!")
            return 1

        factory.emissions_indel = e_indel

    return factory


class AlignSubcommand(Subcommand):
    """
    Align query sequences to a panel of reference sequences.

    If there are multiple query sequences in the FASTA, they will all be aligned to the same reference panel. If you
    want a different reference panel per query contig, use the `multi-align` command.
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

        add_hmm_config_options(subparser)

        output_group = subparser.add_argument_group("Output format")
        output_group.add_argument(
            '-O', '--output-type', default="sam", choices=OUTPUT_TYPES.keys(),
            help="Output type, options include sam or bam output."
        )

        output_group.add_argument(
            '-o', '--output', default="-",
            help="Output filename. Defaults to stdout."
        )

    def __call__(self, ref_panel: Path, query_seqs: Path, output_type: str, output: str, hmm_config=None,
                 *args, **kwargs):  # type: ignore
        logger.info("Configuring HMM...")

        with open_compressed(ref_panel, "rt") as f:
            ref_contigs = list(skbio.io.read(f, "fasta"))

        hmm_params = {}
        if hmm_config is not None:
            config = tomllib.load(hmm_config)
            hmm_params.update(config.get('hmm', {}))

        hmm_params.update({k: v for k, v in kwargs.items() if v is not None})

        factory = configure_hmm_factory(hmm_params)
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


class MultiAlignSubcommand(Subcommand):
    """
    Align multiple query contigs (stored in a single FASTA) to their unique reference panels.

    This command differs from the standard `align` command which aligns all query contigs to the same single reference
    panel, where this command loads up each query contig's unique reference panel given a path template.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'ref_panel_template',
            help="Reference panel path template, to locate the reference panel FASTA for each query contig. Example: "
                 "\"ref_contigs/{query_id}.fa\". Tesserea2 will automatically replace '{query_id}' with the actual "
                 "query contig ID."
        )

        subparser.add_argument(
            'query_seqs', type=Path,
            help="Path to a FASTA file with the query sequences to align. Can be compressed with gzip or bzip2."
        )

        subparser.add_argument(
            '-r', '--range', default=None, required=False, dest='qry_range',
            help="Only align the query contigs within the given range. Format: START:STOP, to align the contigs within "
                 "the zero-based, half open interval [START, STOP). Similar to Python indexing."
        )

        subparser.add_argument(
            '-m', '--max-memory', default=None, required=False,
            help="Maximum amount of memory to allocate. Query sequences that are too long to align within the given "
                 "amount of memory will be skipped and written to a file 'unaligned.fasta.gz'"
        )

        add_hmm_config_options(subparser)

        output_group = subparser.add_argument_group("Output format")
        output_group.add_argument(
            '-O', '--output-type', default="sam", choices=OUTPUT_TYPES.keys(),
            help="Output type, options include sam or bam output."
        )

        output_group.add_argument(
            '-o', '--output-dir', type=Path,
            help="Output directory."
        )

    def __call__(self, ref_panel_template, query_seqs, output_type, output_dir, qry_range=None, hmm_config=None,
                 max_memory=None, *args, **kwargs):
        output_mode = OUTPUT_TYPES[output_type]

        program = {
            'ID': 'tesserae0',
            'PN': 'tesserae',
            'VN': __version__,
            'CL': " ".join(sys.argv)
        }

        hmm_params = {}
        if hmm_config is not None:
            config = tomllib.load(hmm_config)
            hmm_params.update(config.get('hmm', {}))

        hmm_params.update({k: v for k, v in kwargs.items() if v is not None})
        factory = configure_hmm_factory(hmm_params)

        output_dir.mkdir(exist_ok=True, parents=True)

        start = None
        stop = None
        if qry_range:
            parts = qry_range.split(':', maxsplit=1)
            if len(parts) == 1:
                if qry_range.startswith(':'):
                    stop = int(parts[0])
                else:
                    start = int(parts[0])
            else:
                start = int(parts[0])
                stop = int(parts[1])

        if max_memory:
            max_memory = parse_num_suffix(max_memory, bytes=True)

            # Set max_memory to 80% of actual given value to have some buffer
            max_memory = int(round(max_memory * 0.8))

        unaligned_fasta = output_dir / "unaligned.fasta.gz"
        with open_compressed(query_seqs, "rt") as f, open_compressed(unaligned_fasta, "wt") as unaligned_o:
            for i, r in enumerate(skbio.io.read(f, "fasta")):
                if start is not None and i < start:
                    continue

                if stop is not None and i >= stop:
                    continue

                query_id = r.metadata['id']
                query_str = str(r)
                logger.info("Aligning query '%s'...", query_id)

                # Load corresponding reference panel belonging to this query
                ref_panel = Path(ref_panel_template.replace("{query_id}", query_id))
                if not ref_panel.is_file():
                    logger.error("Could not find reference panel FASTA for query '%s'. File '%s' does not exists.",
                                 query_id, ref_panel)
                    continue

                with open_compressed(ref_panel, "rt") as rf:
                    ref_contigs = list(skbio.io.read(rf, "fasta"))

                if not ref_contigs:
                    r.metadata['description'] += " reason=no_refs"
                    skbio.io.write(r, "fasta", into=unaligned_o)
                    logger.info("Skipping query %s, reference panel is empty.", r.metadata['id'])
                    continue

                max_ref_length = max(len(ref) for ref in ref_contigs)

                hmm = factory.build(ref_contigs)

                # Viterbi matrix is float32 (4 bytes) and pointer matrix is two times uint16, so total 8 bytes
                est_mem_usage = len(r) * max_ref_length * hmm.num_states * 8
                if max_memory and est_mem_usage > max_memory:
                    r.metadata['description'] += " reason=too_much_mem"
                    skbio.io.write(r, "fasta", into=unaligned_o)
                    logger.info("Skipping query %s, would require too much memory to align.", r.metadata['id'])
                    continue

                try:
                    result = hmm.align(query_str.encode('ascii'))
                except Exception as e:
                    r.metadata['description'] += " reason=error"
                    skbio.io.write(r, "fasta", into=unaligned_o)

                    logger.error("Error aligning query '%s': %s", query_id, str(e))
                    continue

                output_fname = output_dir / f"{query_id}.{output_type}"

                with SamWriter(ref_contigs, output_fname, output_mode, program) as w:
                    for j, interval in enumerate(result.get_aligned_ref_intervals()):
                        w.write_alignment(r.metadata['id'], query_str, interval, j, result.get_log_likelihood())

        logger.info("Done.")
