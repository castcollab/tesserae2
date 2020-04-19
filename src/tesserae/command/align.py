import argparse
import contextlib
import logging
import math
import os
import pathlib
import sys
import time
from collections import OrderedDict, namedtuple

import pysam

from .. import model

################################################################################

LOGGER = logging.getLogger(__name__)

################################################################################

DEFAULT_OUTPUT_FILE_NAME = "tesserae_alignment.bam"

# Default our quality to zero for now:
DEFAULT_BASE_QUALITY_CHAR = "!"

# Default max PL:
MAX_ALIGNMENT_PL = 60

# Named tuple to store alignment information:
DetailedAlignmentInfo = namedtuple(
    "DetailedAlignmentInfo", ["start_pos", "template_length", "cigar", "qual_pl"]
)

################################################################################


def print_logo():
    """Print the logo to the log."""

    LOGGER.info("=============================================================")
    LOGGER.info("         _____                                 ")
    LOGGER.info("        |_   _|__  ___ ___  ___ _ __ __ _  ___ ")
    LOGGER.info("          | |/ _ \\/ __/ __|/ _ \\ '__/ _` |/ _ \\")
    LOGGER.info("          | |  __/\\__ \\__ \\  __/ | | (_| |  __/")
    LOGGER.info("          |_|\\___||___/___/\\___|_|  \\__,_|\\___|")
    LOGGER.info("                                               ")
    LOGGER.info("")
    LOGGER.info("")
    LOGGER.info("                     +___________+")
    LOGGER.info("                    /:\\         ,:\\")
    LOGGER.info("                   / : \\       , : \\")
    LOGGER.info("                  /  :  \\     ,  :  \\")
    LOGGER.info("                 /   :   +-----------+")
    LOGGER.info("                +....:../:...+   :  /|")
    LOGGER.info("                |\\   +./.:...`...+ / |")
    LOGGER.info("                | \\ ,`/  :   :` ,`/  |")
    LOGGER.info("                |  \\ /`. :   : ` /`  |")
    LOGGER.info("                | , +-----------+  ` |")
    LOGGER.info("                |,  |   `+...:,.|...`+")
    LOGGER.info("                +...|...,'...+  |   /")
    LOGGER.info("                 \\  |  ,     `  |  /")
    LOGGER.info("                  \\ | ,       ` | /")
    LOGGER.info("                   \\|,         `|/ ")
    LOGGER.info("                    +___________+")
    LOGGER.info("")
    LOGGER.info("=============================================================")
    LOGGER.info("")


def compute_detailed_alignment_info(
    query_alignment_string, target_alignment_string, target_length
):
    """Compute detailed alignment information from the given information.

    Alignment details are based off the differences between the alignment
    strings.
    This method returns a tuple containing:
    - The Start Position in the reference of the alignment.
    - The Template Length of this alignment.
    - The Cigar representing this alignment.
    - The Phred-Scaled quality score of this alignment.
    Where:
    - The Start Position is the 1-based, inclusive position in the reference
    at which this alignment begins.
    - The Template Length is the number of bases accounted by this alignment
    with respect to the reference.
    - The Cigar is a list of tuples: (CIGAR_ELEMENT, COUNT) where each
    CIGAR_ELEMENT is defined in pysam.
    - The Phred-Scaled quality score is defined by the following formula:
    -10 log_10((# mismatches + # insertions + # deletions)/target_length)
    """

    start_index = get_start_index_from_alignment_start_string(target_alignment_string)

    # Now that we know where in the reference this target begins, we can start
    # to loop through both alignment strings at the same time.
    # In this loop we will:
    #     construct a cigar string
    #     determine counts for alignment quality score
    #     determine template length

    num_errors = 0

    cigar = []
    current_cigar_element = None
    current_cigar_element_count = 0

    for query_base, target_base in zip(
        query_alignment_string[start_index:], target_alignment_string[start_index:]
    ):

        # The Tesserae2 / Mosaic Alignment algorithm can only produce "-" or
        # <BASE> for any position (other than blanks / spaces).  Therefore we
        # only have to check the following 4 cases:
        if query_base == "-":
            # We have an insertion relative to the reference:
            num_errors += 1
            cigar_element = pysam.CINS

        elif query_base == target_base:
            # Bases match:
            # We use CMATCH here because that cigar operator accounts for
            # BOTH matches and mismatches.
            cigar_element = pysam.CMATCH

        elif target_base == "-":
            # We have a deletion relative to the reference:
            num_errors += 1
            cigar_element = pysam.CDEL

        else:
            # We have a mismatch relative to the reference:
            num_errors += 1
            # We use CMATCH here because that cigar operator accounts for
            # BOTH matches and mismatches.
            cigar_element = pysam.CMATCH

        # Accumulate our cigar elements:
        if cigar_element != current_cigar_element:
            if current_cigar_element is not None:
                cigar.append((current_cigar_element, current_cigar_element_count))

            current_cigar_element = cigar_element
            current_cigar_element_count = 1

        else:
            current_cigar_element_count += 1

    # Add the last remaining cigar element to our list:
    cigar.append((current_cigar_element, current_cigar_element_count))

    # Our template length is the number of bases accounted by this alignment
    # with respect to the reference:
    template_length = len(target_alignment_string) - start_index

    # Compute PL score:
    if num_errors == 0:
        qual_pl = MAX_ALIGNMENT_PL
    else:
        qual_pl = -10 * math.log10(num_errors / target_length)
    if qual_pl < 0:
        # quality cannot take a negative value
        qual_pl = 0

    # Return our detailed info.
    return DetailedAlignmentInfo(start_index, template_length, cigar, qual_pl)


def get_start_index_from_alignment_start_string(target_alignment_string):
    """Get the alignment start index from the given alignment string.

    We know that the alignment strings will always start with spaces until
    the region that actually aligns to the reference, so we count the number
    of spaces in the target_alignment_string to get our starting alignment
    position.
    """

    start_index = 0
    while start_index < len(target_alignment_string):
        if target_alignment_string[start_index] != " ":
            break
        start_index += 1
    return start_index


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
    with pysam.FastxFile(file_path) as file_handle:
        for entry in file_handle:
            _read_to_sequence_dict[entry.name] = entry.sequence
            _template_to_read_name_dict[f"template{t_num}"] = entry.name

            t_num += 1

    return _read_to_sequence_dict, _template_to_read_name_dict


def dump_results_to_log(results_tuple_list):
    """Dump the results_tuple_list to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Results:")
        for result in results_tuple_list:
            LOGGER.debug("    %s", str(result))


def write_results(query, list_of_result_tuples, bamfile=None):
    """Write the alignment results in the given `list_of_result_tuples` to the
    stdout and to the given `bamfile` (if present).

    Assumes that list_of_result_tuples does NOT contain information from the
    Query sequence itself.

    NOTE: If a file already exists a the given `bamfile` location, the file will be
    overwritten.

    If the path to the given bamfile does not exist, will warn the user and
    output the alignment information in a file of the given name in the current
    working directory.

    Each tuple in the given `list_of_result_tuples` should have the following
    fields:

        - target
        - alignment start_pos
        - template_length
        - cigar
        - alignment_qual_pl
        - target_start_index
        - target_end_index

    """

    bamfile = prepare_output_file(bamfile)

    # Create a header for the output bam file wiith the query read as the
    # reference.
    # (The query should be the first entry in the tuple list.)
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": len(query), "SN": query.name}],
    }

    # Output the sam results to stdout and to a file if the user specified it.
    sam_stdout = pysam.AlignmentFile("-", "w", header=header)
    with contextlib.ExitStack() as stack:
        stack.enter_context(contextlib.closing(sam_stdout))

        bam_output_file = None
        if bamfile:
            bam_output_file = pysam.AlignmentFile(bamfile, "wb", header=header)
            stack.enter_context(contextlib.closing(bam_output_file))

        for result in list_of_result_tuples:

            # Unpack our tuple for ease-of-use:
            (
                target,
                ref_start_pos,
                template_length,
                cigar,
                alignment_qual_pl,
                target_start_index,
                target_end_index,
            ) = result

            # Create our segment:
            aligned_segment = pysam.AlignedSegment()

            # Populate the name:
            aligned_segment.query_name = (
                f"{target.name}_{target_start_index}_{target_end_index}"
            )

            # We have only one reference, so use it:
            aligned_segment.reference_id = 0

            # Set our sequence being aligned:
            # NOTE: The nomenclature here is overloaded. Query in this case
            # means the sequence being aligned
            #       to the reference.
            aligned_segment.query_sequence = target.sequence[
                target_start_index : target_end_index + 1
            ]

            aligned_segment.reference_start = ref_start_pos
            aligned_segment.template_length = template_length
            aligned_segment.cigar = cigar
            aligned_segment.mapping_quality = alignment_qual_pl

            aligned_segment.query_qualities = pysam.qualitystring_to_array(
                DEFAULT_BASE_QUALITY_CHAR * len(aligned_segment.query_sequence)
            )

            if bam_output_file:
                bam_output_file.write(aligned_segment)

            sam_stdout.write(aligned_segment)


def prepare_output_file(bamfile):
    """Prepare the given bamfile path to be ready to accept data.

    The input `bamfile` can be a file or a directory.

    If the given `bamfile` is a file and does not exist, the path to the file is
    created.
    If the given `bamfile` already exists a warning is logged.
    If the given `bamfile` is a path, the DEFAULT_OUTPUT_FILE_NAME is appended to
    this path and then used as the bamfile location.

    The final bamfile path is returned as a string object.
    """
    # Do some work to ensure that the given file will be honored.
    # Since we do this at the end of execution, we should be liberal with our
    # output file creation to preserve the work that has happened already:

    if bamfile is not None:
        out_path = pathlib.Path(bamfile)
        if out_path.is_file():
            # File already exists:
            LOGGER.warning("Overwriting existing file: %s", bamfile)
        elif out_path.is_dir():
            # File path is actually a directory:
            LOGGER.warning("Given output file path is a directory: %s", bamfile)

            # Create a proper location for an output file:
            bamfile = bamfile + os.path.sep + DEFAULT_OUTPUT_FILE_NAME

            LOGGER.warning("Creating output file: %s", bamfile)
        elif out_path.parent.absolute().is_dir():
            # File path does not exist, but is in an existing directory:

            LOGGER.warning("Creating output file: %s", bamfile)

        elif not out_path.parent.absolute().exists():
            # File path parent does not exist.

            # Make containing directory:
            LOGGER.info(
                "Creating containing directory for output file: %s",
                out_path.parent.absolute(),
            )
            os.makedirs(out_path.parent.absolute(), exist_ok=True)

            LOGGER.warning("Creating output file: %s", bamfile)

    return bamfile


################################################################################


def align(args):
    """Main CLI call for the 'align' sub-command."""

    aligner = model.Tesserae()
    LOGGER.info("Aligning with Tesserae...")
    start_time = time.time()
    target_alignment_results = aligner.align_from_fastx(args.reads, args.segments)
    end_time = time.time()
    LOGGER.info("Alignment took %fs", end_time - start_time)
    dump_results_to_log(target_alignment_results)

    LOGGER.info("Cleaning results and computing alignment details...")
    query_alignment_string = target_alignment_results[0].alignment_string
    clean_results = []
    for result in target_alignment_results[1:]:
        # Unpack our results for ease of use:

        # Get our detailed alignment info:
        target = aligner.get_target(result.seq_name)
        detailed_alignment_info = compute_detailed_alignment_info(
            query_alignment_string, result.alignment_string, len(target)
        )

        # Create our clean results for this alignment:
        clean_results.append(
            (
                target,
                detailed_alignment_info.start_pos,
                detailed_alignment_info.template_length,
                detailed_alignment_info.cigar,
                detailed_alignment_info.qual_pl,
                result.target_start_index,
                result.target_end_index,
            )
        )

    dump_results_to_log(clean_results)

    # Get our query information:
    query = aligner.query

    LOGGER.info("Writing results...")
    write_results(query, clean_results, args.bamout)

    LOGGER.debug("Dumping raw tesserae object:")
    LOGGER.debug(aligner)


################################################################################


def main(raw_args):

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        prog="tesserae align",
        description="Reads two input fasta files - a reads file and a known "
        "segments file.  The known segments are aligned to each "
        "read in the reads file and the resulting alignments are "
        "output in SAM file format on stdout (and optionally to a file).",
        usage="align query and target sequences",
    )

    align_required_args = parser.add_argument_group("Required Arguments")
    align_required_args.add_argument(
        "-r", "--reads", help="Reads FASTA filename.", required=True
    )
    align_required_args.add_argument(
        "-s", "--segments", help="Known segments FASTA filename.", required=True
    )

    parser.add_argument(
        "-o",
        "--bamout",
        help="Output BAM file in which to store alignment results.",
        required=False,
    )

    # ---------------------------------------

    # Parse args
    args = parser.parse_args(args=raw_args)

    # Print logo:
    print_logo()

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(sys.argv))
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our sub-command:
    align(args)

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)
