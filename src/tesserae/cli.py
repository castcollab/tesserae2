import time
import pysam
import argparse
import math

import os
import sys
import pathlib

from collections import OrderedDict

from . import tesserae
from . import log

################################################################################

logger = log.get_logger(__name__)

################################################################################

DEFAULT_OUTPUT_FILE_NAME = "integration_test_expected_output.bam"

# Default our quality to zero for now:
DEFAULT_BASE_QUALITY_CHAR = "!"

# Default max PL:
MAX_ALIGNMENT_PL = 60

################################################################################


def print_logo():
    """Prints the logo to the log."""

    logger.info("=============================================================")
    logger.info("         _____                                 ")
    logger.info("        |_   _|__  ___ ___  ___ _ __ __ _  ___ ")
    logger.info("          | |/ _ \\/ __/ __|/ _ \\ '__/ _` |/ _ \\")
    logger.info("          | |  __/\\__ \\__ \\  __/ | | (_| |  __/")
    logger.info("          |_|\\___||___/___/\\___|_|  \\__,_|\\___|")
    logger.info("                                               ")
    logger.info("")
    logger.info("")
    logger.info("                     +___________+")
    logger.info("                    /:\\         ,:\\")
    logger.info("                   / : \\       , : \\")
    logger.info("                  /  :  \\     ,  :  \\")
    logger.info("                 /   :   +-----------+")
    logger.info("                +....:../:...+   :  /|")
    logger.info("                |\\   +./.:...`...+ / |")
    logger.info("                | \\ ,`/  :   :` ,`/  |")
    logger.info("                |  \\ /`. :   : ` /`  |")
    logger.info("                | , +-----------+  ` |")
    logger.info("                |,  |   `+...:,.|...`+")
    logger.info("                +...|...,'...+  |   /")
    logger.info("                 \\  |  ,     `  |  /")
    logger.info("                  \\ | ,       ` | /")
    logger.info("                   \\|,         `|/ ")
    logger.info("                    +___________+")
    logger.info("")
    logger.info("=============================================================")
    logger.info("")


def compute_detailed_alignment_info(
    query_alignment_string, target_alignment_string, target_length
):
    """Computes detailed alignment information from the given information.

    Alignment details are based off the differences between the alignment strings.
    This method returns a tuple containing:
        - The Start Position in the reference of the alignment.
        - The Template Length of this alignment.
        - The Cigar representing this alignment.
        - The Phred-Scaled quality score of this alignment.
    Where:
        - The Start Position is the 1-based, inclusive position in the reference at which this alignment begins.
        - The Template Length is the number of bases accounted by this alignment with respect to the reference.
        - The Cigar is a list of tuples: (CIGAR_ELEMENT, COUNT) where each CIGAR_ELEMENT is defined in pysam.
        - The Phred-Scaled quality score is defined by the following formula:
            -10 log_10((# mismatches + # insertions + # deletions)/target_length)
    """

    # We know that the alignment strings will always start with spaces until the region that actually aligns to the
    # reference, so we count the number of spaces in the target_alignment_string to get our starting alignment
    # position:
    start_index = 0
    while start_index < len(target_alignment_string):
        if target_alignment_string[start_index] != " ":
            break
        start_index += 1

    # Now that we know where in the reference this target begins, we can start to loop through both alignment
    # strings at the same time.
    # In this loop we will:
    #     construct a cigar string
    #     determine counts for alignment quality score
    #     determine template length

    num_mismatches = 0
    num_insertions = 0
    num_deletions = 0
    num_matches = 0

    cigar = []
    current_cigar_element = None
    current_cigar_element_count = 0

    i = start_index
    while (i < len(query_alignment_string)) and (i < len(target_alignment_string)):

        # The Tesserae2 / Mosaic Alignment algorithm can only produce "-" or <BASE> for any position
        # (other than blanks / spaces).  Therefore we only have to check the following 4 cases:
        if query_alignment_string[i] == "-":
            # We have an insertion relative to the reference:
            num_insertions += 1
            cigar_element = pysam.CINS

        elif query_alignment_string[i] == target_alignment_string[i]:
            # Bases match:
            num_matches += 1
            # We use CMATCH here because that cigar operator accounts for BOTH matches and mismatches.
            cigar_element = pysam.CMATCH

        elif target_alignment_string[i] == "-":
            # We have a deletion relative to the reference:
            num_deletions += 1
            cigar_element = pysam.CDEL

        else:
            # We have a mismatch relative to the reference:
            num_mismatches += 1
            # We use CMATCH here because that cigar operator accounts for BOTH matches and mismatches.
            cigar_element = pysam.CMATCH

        # Accumulate our cigar elements:
        if cigar_element != current_cigar_element:
            if current_cigar_element is not None:
                cigar.append((current_cigar_element, current_cigar_element_count))

            current_cigar_element = cigar_element
            current_cigar_element_count = 1

        else:
            current_cigar_element_count += 1

        # Advance to the next base:
        i += 1

    # Make sure we account for any remaining bases that have not been accounted for:
    if (i < len(query_alignment_string)) or (i < len(target_alignment_string)):
        if i >= len(query_alignment_string):
            # Having extra bases at the end of the query is OK
            pass
        if i >= len(target_alignment_string):
            # Any extra bases at the end of the target alignment have been deleted relative to the query
            # and must be included in our cigar:
            if current_cigar_element != pysam.CDEL:
                # Add in the remaining cigar elements:
                cigar.append((current_cigar_element, current_cigar_element_count))

                # Set the remaining bases as an insertion so we can add them below:
                current_cigar_element = pysam.CDEL
                current_cigar_element_count = i - len(target_alignment_string)
            else:
                current_cigar_element_count += i - len(target_alignment_string)

    # Add the last remaining cigar element to our list:
    cigar.append((current_cigar_element, current_cigar_element_count))

    # Our template length is the number of bases accounted by this alignment with respect to the reference:
    template_length = len(target_alignment_string) - start_index

    # Compute PL score:
    num_errors = num_mismatches + num_insertions + num_deletions
    if num_errors == 0:
        quality_pl = MAX_ALIGNMENT_PL
    else:
        quality_pl = -10 * math.log10(num_errors / target_length)

    # Return our detailed info.
    # NOTE: The start position in BAM / SAM files is 1-based.
    return start_index + 1, template_length, cigar, quality_pl


def ingest_fastx_file(file_path):
    """Ingest the contents of a FASTA/FASTQ file and return two dictionaries:
        1: Mapping from read name to read sequence
        2: Mapping from template name to read name

    The `template name` is simply the word 'template' with the 
    order in which a given read occurs in the file 
    (e.g. 'template0' or 'template10').
    """

    t_num = 0

    _read_to_sequence_dict = OrderedDict()
    _template_to_read_name_dict = dict()
    with pysam.FastxFile(file_path) as fh:
        for _entry in fh:
            _read_to_sequence_dict[_entry.name] = _entry.sequence
            _template_to_read_name_dict["template" + str(t_num)] = _entry.name

            t_num += 1

    return _read_to_sequence_dict, _template_to_read_name_dict


def dump_seq_map(seq_map, name="Reads"):
    """Dumps the given sequence map to the log as a DEBUG message."""
    if logger.isEnabledFor(log.DEBUG):
        logger.debug("%s:", name)
        for e in seq_map.items():
            logger.debug("    %s -> %s", e[0], e[1])


def dump_results(results_tuple_list):
    """Dumps the results tuple to the log as a DEBUG message."""
    if logger.isEnabledFor(log.DEBUG):
        logger.debug("Results:")
        for r in results_tuple_list:
            logger.debug("    %s", str(r))


def write_results_to_bam_file(
    bamfile, query_name, query_sequence, list_of_result_tuples
):
    """Writes the alignment results in the given `list_of_result_tuples` to the given `bamfile`.
    Assumes that list_of_result_tuples does NOT contain information from the Query sequence itself.

    NOTE: If a file already exists a the given output location, the file will be overwritten.

    If the path to the given bamfile does not exist, will warn the user and output the alignment information
    in a file of the given name in the current working directory.

    Each tuple in the given `list_of_result_tuples` should have the following fields:

        - target_name
        - target_sequence
        - alignment start_pos
        - template_length
        - cigar
        - alignment_qual_pl
        - target_start_index
        - target_end_index

    """

    # Do some work to ensure that the given file will be honored.
    # Since we do this at the end of execution, we should be liberal with our output file creation
    # to preserve the work that has happened already:
    out_path = pathlib.Path(bamfile)

    if out_path.is_file():
        # File already exists:
        logger.warning("Overwriting existing file: " + bamfile)
    elif out_path.is_dir():
        # File path is actually a directory:
        logger.warning("Given output file path is a directory: " + bamfile)

        # Create a proper location for an output file:
        bamfile = bamfile + os.path.sep + DEFAULT_OUTPUT_FILE_NAME

        logger.warning("Creating output file: " + bamfile)
    elif out_path.parent.absolute().is_dir():
        # File path does not exist, but is in an existing directory:

        logger.warning("Creating output file: " + bamfile)

    elif not out_path.parent.absolute().exists():
        # File path parent does not exist.

        # Make containing directory:
        logger.info(
            "Creating containing directory for output file: "
            + str(out_path.parent.absolute())
        )
        os.makedirs(out_path.parent.absolute(), exist_ok=True)

        logger.warning("Creating output file: " + bamfile)

    # Create a header for the output bam file wiith the query read as the reference.
    # (The query should be the first entry in the tuple list.)
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": len(query_sequence), "SN": query_name}],
    }

    with pysam.AlignmentFile(bamfile, "wb", header=header) as bam_out:
        for result in list_of_result_tuples:

            # Unpack our tuple for ease-of-use:
            (
                target_name,
                target_sequence,
                ref_start_pos,
                template_length,
                cigar,
                alignment_qual_pl,
                target_start_index,
                target_end_index,
            ) = result

            # Create our segment:
            a = pysam.AlignedSegment()

            # Populate the name:
            a.query_name = (
                target_name
                + "_"
                + str(target_start_index)
                + "_"
                + str(target_end_index)
            )

            # We have only one reference, so use it:
            a.reference_id = 0

            # Set our sequence being aligned:
            # NOTE: The nomenclature here is overloaded. Query in this case means the sequence being aligned
            #       to the reference.
            a.query_sequence = target_sequence[
                target_start_index : target_end_index + 1
            ]

            a.reference_start = ref_start_pos
            a.template_length = template_length
            a.cigar = cigar
            a.mapping_quality = alignment_qual_pl

            a.query_qualities = pysam.qualitystring_to_array(
                DEFAULT_BASE_QUALITY_CHAR * len(a.query_sequence)
            )

            bam_out.write(a)


################################################################################


def align(_args):
    """Main CLI call for the 'align' sub-command."""

    t = tesserae.Tesserae()
    logger.info("Aligning with Tesserae...")
    start_time = time.time()
    results = t.align_from_fastx(_args.reads, _args.segments)
    end_time = time.time()
    logger.info("Alignment took %fs", end_time - start_time)
    dump_results(results)

    logger.info("Cleaning results and computing alignment details...")
    query_alignment_string = results[0][1]
    clean_results = []
    for i in range(1, len(results)):
        # Unpack our results for ease of use:
        target_name, aligned_sequence, target_start_index, target_end_index = results[i]

        # Get our detailed alignment info:
        target_sequence = t.targets[target_name]
        target_length = len(target_sequence)
        (
            ref_start_pos,
            template_length,
            cigar,
            alignment_qual_pl,
        ) = compute_detailed_alignment_info(
            query_alignment_string, aligned_sequence, target_length
        )

        # Create our clean results for this alignment:
        clean_results.append(
            (
                target_name,
                target_sequence,
                ref_start_pos,
                template_length,
                cigar,
                alignment_qual_pl,
                target_start_index,
                target_end_index,
            )
        )

    # Sort our results (by sequence name):
    sorted(clean_results)

    # Dump our raw results to the log:
    dump_results(clean_results)

    # Get our query information:
    query_name = list(t.query.keys())[0]
    query_sequence = list(t.query.values())[0]

    # Write our results to the output BAM file:
    logger.info("Writing results...")
    write_results_to_bam_file(_args.bamout, query_name, query_sequence, clean_results)

    # Some debug logging:
    logger.debug("Dumping raw tesserae object:")
    logger.debug(t)


################################################################################


def main():
    """Main CLI entry point for Tesserae"""

    # Get our start time:
    overall_start = time.time()

    # Setup main parser:
    parser = argparse.ArgumentParser(
        prog="tesserae",
        description="Graph-based mosaic read alignment and exploration algorithms.",
        epilog="",
    )
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )

    # Prepare for sub-parsers:
    subparsers = parser.add_subparsers(
        title="sub-commands", dest="sub-command", help=""
    )
    subparsers.required = True

    # ---------------------------------------
    # Setup parser for 'align':
    align_parser = subparsers.add_parser(
        "align",
        description="Reads two input fasta files - a reads file and a known segments file.  The known segments are "
        + "aligned to each read in the reads file and the resulting alignments are output as a sam file.  "
        + "Unaccounted for regions in the reads that are not accounted for in the alignments of the known "
        + "segments are also identified and output in the resulting file.",
        help="align query and target sequences",
        epilog="",
    )

    align_required_args = align_parser.add_argument_group("Required Arguments")
    align_required_args.add_argument(
        "-r", "--reads", help="Reads FASTA filename.", required=True
    )
    align_required_args.add_argument(
        "-s", "--segments", help="Known segments FASTA filename.", required=True
    )

    align_parser.add_argument(
        "-o",
        "--bamout",
        help="Output BAM file in which to store alignment results.",
        required=False,
        default=DEFAULT_OUTPUT_FILE_NAME,
    )

    # Set our subcommand to point to the `align` function:
    align_parser.set_defaults(func=align)

    # ---------------------------------------

    # Parse args
    args = parser.parse_args()

    # Set logging level:
    if args.verbose:
        log.set_level(log.DEBUG)

    # Print logo:
    print_logo()

    # Log our command-line so we can have it in the log file:
    logger.info("Invoked by: %s", " ".join(sys.argv))

    # Call our sub-command:
    args.func(args)

    overall_end = time.time()
    logger.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == "__main__":
    main()
