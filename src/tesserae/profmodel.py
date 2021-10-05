import logging
import re
import math

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple, Optional, Any

import pysam
from pomegranate import *

from .nucleotide_sequence import NucleotideSequence

LOGGER = logging.getLogger(__name__)

# default transition probabilities
DEFAULT_DEL = 0.025
DEFAULT_EPS = 0.75
DEFAULT_REC = 0.0001
DEFAULT_TERM = 0.001

MAX_ALIGNMENT_PL = 60

# by default do no subsampling (every site is considered for recombination)
DEFAULT_SUBSAMPLING = 1

# Nucleotide to index conversion
# fmt: off
CONVERT = {
    'G': 4, 'g': 4,
    'A': 3, 'a': 3,
    'C': 2, 'c': 2,
    'T': 1, 't': 1, 'U': 1, 'u': 1,
    'N': 0, 'n': 0,
}
# fmt: on

@dataclass
class DetailedAlignmentInfo:
    """Stores alignment information."""

    start_pos: int
    template_length: int
    cigar: List[Tuple[Optional[Any], int]]
    qual_pl: float


@dataclass
class TesseraeAlignmentResult:
    """
    Representation of the results of a Tesserae alignment
    """

    seq_name: str
    alignment_string: str
    target_start_index: int
    target_end_index: int

@dataclass
class CleanResult:
    """Clean alignment Result summarizing the alignment of the target to a particular sequence with summarizing information."""

    target: NucleotideSequence
    ref_start_pos: int
    template_length: int
    cigar: Tuple
    alignment_qual_pl: float
    target_start_index: int
    target_end_index: int

    @classmethod
    def from_target__detailed_alignment_info_and_target_alignment_result(
            cls, target, detailed_alignment_info, target_alignment_result
    ):
        return cls(
            target,
            detailed_alignment_info.start_pos,
            detailed_alignment_info.template_length,
            detailed_alignment_info.cigar,
            detailed_alignment_info.qual_pl,
            target_alignment_result.target_start_index,
            target_alignment_result.target_end_index,
        )

    def get_query_name(self) -> str:
        """Build the name of the query sequence"""
        return "_".join(
            str(v)
            for v in [self.target.name, self.target_start_index, self.target_end_index]
        )

@dataclass
class ProfalignAlignmentResultContainer:
    """
    Class for storing the necessary Tesserae results
    """

    aligner: Any
    pomegranate_path: List[Any]
    alignment_result_list: List[TesseraeAlignmentResult]
    clean_result_summary: List[CleanResult]

    logp: float

    def get_pomegranate_path(self) -> List[Any]:
        return self.pomegranate_path.copy()

    def get_blat_results(self) -> List[TesseraeAlignmentResult]:
        return self.alignment_result_list.copy()

    def get_clean_alignment_results(self) -> List[CleanResult]:
        query_alignment_string = self.alignment_result_list[0].alignment_string
        self.clean_result_summary = clean_results = []

        for result in self.alignment_result_list[1:]:
            # Unpack our results for ease of use:

            if not result.seq_name.startswith("flank") and not result.seq_name.startswith("recombination"):
                # Get our detailed alignment info:
                target = self.aligner.get_target(result.seq_name)
                detailed_alignment_info = self.compute_detailed_alignment_info(
                    query_alignment_string, result.alignment_string, len(target)
                )
            else:
                target = NucleotideSequence(result.seq_name, result.alignment_string)
                detailed_alignment_info = self.compute_flanking_alignment_info(
                    query_alignment_string, result.alignment_string
                )
                if detailed_alignment_info.template_length is 0:
                    continue

            # Create our clean results for this alignment:
            clean_results.append(
                CleanResult.from_target__detailed_alignment_info_and_target_alignment_result(
                    target, detailed_alignment_info, result
                )
            )
        return self.clean_result_summary

    def compute_detailed_alignment_info(
            cls, query_alignment_string, target_alignment_string, target_length
    ) -> DetailedAlignmentInfo:
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

        start_index = cls.get_start_index_from_alignment_start_string(target_alignment_string)

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
        qual_pl: float
        if num_errors == 0:
            qual_pl = MAX_ALIGNMENT_PL
        else:
            qual_pl = -10 * math.log10(num_errors / target_length)
        if qual_pl < 0:
            # quality cannot take a negative value
            qual_pl = 0

        # Return our detailed info.
        return DetailedAlignmentInfo(start_index, template_length, cigar, qual_pl)

    def compute_flanking_alignment_info(cls,
                                        query_alignment_string, target_alignment_string
                                        ) -> DetailedAlignmentInfo:
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

        start_index = cls.get_start_index_from_alignment_start_string(target_alignment_string)

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

        for query_xbase, target_base in zip(
                query_alignment_string[start_index:], target_alignment_string[start_index:]
        ):
            cigar_element = pysam.CSOFT_CLIP

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
        qual_pl = 0.0

        # Return our detailed info.
        return DetailedAlignmentInfo(start_index, template_length, cigar, qual_pl)

    def get_start_index_from_alignment_start_string(self, target_alignment_string) -> int:
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


def repeat(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]


def ingest_fastx_file(fastx_file_name: str) -> List[NucleotideSequence]:
    """Ingest the given fastx (fasta or fastq) file and returns a
    list of Sequence objects."""

    LOGGER.info("Ingesting reads from %s ...", fastx_file_name)

    sequence_list = []
    with pysam.FastxFile(fastx_file_name) as fh:
        for entry in fh:
            sequence_list.append(NucleotideSequence(entry.name, entry.sequence))

    LOGGER.info("Ingested %d reads.", len(sequence_list))

    return sequence_list


class Tesserae2:
    """Object to perform profile HMM-based mosaic alignments between nucleic acid
    sequences.

    All parameters except for query and sources (at the moment) need to be
    specified at initialization time. In other words, all attributes should
    be considered read-only! Anything else may be madness.
    """

    def __init__(
            self, threads=1, subsample_model_recombination=DEFAULT_SUBSAMPLING, use_flanking_model=True
    ):
        self.logp = 0.0
        if sys.version_info < (3, 8) and threads > 1:
            raise NotImplementedError(
                "Threading is only supported with python3.8 and above!"
            )
        self.threads: int = threads
        self.query: NucleotideSequence
        self.subsample_model_recombination: int = subsample_model_recombination
        self.sources: Sequence[NucleotideSequence]
        self.use_flanking_model: bool = use_flanking_model

        self._source_name_index_dict: Dict[str, int] = {}

    def _query_and_source_dependent_setup(self, query, sources, sourceKeys):
        """Perform all initialization that depends on the query or sources."""
        self.query = query
        self.sources = sources
        self.sourceKeys = sourceKeys
        self.model = self._make_full_model()

    def _make_global_alignment_model(self, source, keyArray):
        name = source.name
        pos1extra = "*" if keyArray[1] else ""

        # enforce that the begining and the end of the arrays are always reachable
        keyArray[0] = keyArray[len(keyArray)-1] = True

        model = HiddenMarkovModel(name=name)
        s = {}

        # add states
        i0 = State(
            DiscreteDistribution({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}),
            name=f"{name}:I*0",
        )

        model.add_state(i0)

        s[i0.name] = i0

        for c in range(len(source)):
            extra = "*" if c + 1 >= len(keyArray) or keyArray[c + 1] else ""
            dc = State(None, name=f"{name}:D{extra}{c + 1}")

            # TODO: replace this with the nucleotide components of emission matrix from
            #  model.py (line 148)
            mc = State(
                DiscreteDistribution(
                    {
                        "A": 0.94 if source.sequence[c] == "A" else 0.02,
                        "C": 0.94 if source.sequence[c] == "C" else 0.02,
                        "G": 0.94 if source.sequence[c] == "G" else 0.02,
                        "T": 0.94 if source.sequence[c] == "T" else 0.02,
                    }
                ),
                name=f"{name}:M{extra}{c + 1}",
            )

            ic = State(
                DiscreteDistribution({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}),
                name=f"{name}:I{extra}{c + 1}",
            )

            model.add_states([mc, ic, dc])

            s[dc.name] = dc
            s[mc.name] = mc
            s[ic.name] = ic

        # add transitions
        model.add_transition(model.start, s[f"{name}:I*0"], 0.05)
        model.add_transition(model.start, s[f"{name}:D{pos1extra}1"], 0.05)
        model.add_transition(model.start, s[f"{name}:M{pos1extra}1"], 0.90)

        model.add_transition(s[f"{name}:I*0"], s[f"{name}:I*0"], 0.70)
        model.add_transition(s[f"{name}:I*0"], s[f"{name}:D{pos1extra}1"], 0.15)
        model.add_transition(s[f"{name}:I*0"], s[f"{name}:M{pos1extra}1"], 0.15)

        for c in range(1, len(source)):
            extra = "*" if keyArray[c] else ""
            extraPlus = "*" if c + 1 >= len(keyArray) or keyArray[c + 1] else ""
            num = extra + str(c)
            numPlus = extraPlus + str(c + 1)
            model.add_transition(s[f"{name}:D{num}"], s[f"{name}:D{numPlus}"], 0.15)
            model.add_transition(s[f"{name}:D{num}"], s[f"{name}:I{num}"], 0.70)
            model.add_transition(s[f"{name}:D{num}"], s[f"{name}:M{numPlus}"], 0.15)

            model.add_transition(s[f"{name}:I{num}"], s[f"{name}:D{numPlus}"], 0.15)
            model.add_transition(s[f"{name}:I{num}"], s[f"{name}:I{num}"], 0.15)
            model.add_transition(s[f"{name}:I{num}"], s[f"{name}:M{numPlus}"], 0.70)

            model.add_transition(s[f"{name}:M{num}"], s[f"{name}:D{numPlus}"], 0.05)
            model.add_transition(s[f"{name}:M{num}"], s[f"{name}:I{num}"], 0.05)
            model.add_transition(s[f"{name}:M{num}"], s[f"{name}:M{numPlus}"], 0.90)

        model.add_transition(s[f"{name}:D*{len(source)}"], s[f"{name}:I*{len(source)}"], 0.70)

        model.add_transition(s[f"{name}:I*{len(source)}"], s[f"{name}:I*{len(source)}"], 0.15)

        model.add_transition(s[f"{name}:M*{len(source)}"], s[f"{name}:I*{len(source)}"], 0.90)

        ## add the recombination silent state
        model.add_state(State(None, name="recombination:RC"))

        model.bake(merge="None")

        return model

    def subsample_model_recombination_sites(self, source, subsample_interval):
        array = [idx < 5 or idx + 5 > len(source) or (idx % subsample_interval) == 0 for idx in range(len(source))]
        return array

    def _make_flanking_model(self, name):
        model = HiddenMarkovModel(name=name)

        # add states
        fi = State(
            DiscreteDistribution({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}),
            name=f"{name}:FI",
        )
        fd = State(None, name=f"{name}:FD")

        model.add_states([fi, fd])

        # add transitions
        model.add_transition(model.start, fd, 0.5)
        model.add_transition(model.start, fi, 0.5)

        model.add_transition(fi, fi, 0.8)
        model.add_transition(fi, fd, 0.2)

        model.add_transition(fd, model.end, 1.0)

        model.bake(merge="None")

        return model

    def _make_full_model(self):

        if self.use_flanking_model:
            full_model = self._make_flanking_model("flankleft")
            full_model.add_model(self._make_flanking_model("flankright"))
        else:
            full_model = HiddenMarkovModel("basemodel")

        for i in range(len(self.sources)):
            s = self.sources[i]
            if self.sourceKeys is None:
                source_key_to_use = self.subsample_model_recombination_sites(s, self.subsample_model_recombination)
            else:
                source_key_to_use = self.sourceKeys[i]

            assert len(source_key_to_use) == len(s.sequence), "Size mismatch for provided source recomination site array " \
                                                     "expected " + len(s.sequence) + " but found " + len(source_key_to_use) + " for source:" + s.name

            full_model.add_model(self._make_global_alignment_model(s, source_key_to_use))

        full_model.bake(merge="None")

        # dictionary of model states for easy lookup later
        states = {}
        for state in full_model.states:
            model, name = re.split("[:-]", state.name)

            if model not in states:
                states[model] = {}

            if name == "RC":
                recomb = state

            elif name != "start" and name != "end" and not name.startswith("D"):
                if "*" in name or "FD" in name or "FI" in name:
                    states[model][name] = state

        # link M and I states between different models to each other (i.e. allow for recombination)
        xlink = 0
        for s1 in self.sources:
            for state1 in states[s1.name]:
                full_model.add_transition(states[s1.name][state1], recomb, DEFAULT_REC)
                if "F" not in state1:
                    full_model.add_transition(recomb, states[s1.name][state1], 1)  # only penalize recombinations once
                    xlink += 1
                xlink += 1

        # link flanking models to match states
        if self.use_flanking_model:
            for model in states:
                if not model.startswith("flank"):
                    for s in [x for x in states[model] if x.startswith("M")]:
                        full_model.add_transition(states['flankleft']['FD'], states[model][s], 0.001)
                        full_model.add_transition(states[model][s], states['flankright']['FI'], 0.001)

        # link model start and end to all match and insert states
        for model in states:
            if not model.startswith("flank"):
                for s in [x for x in states[model] if x.startswith("M") or x.startswith("I")]:
                    full_model.add_transition(full_model.start, states[model][s], 0.001)
                    full_model.add_transition(states[model][s], full_model.end, 0.001)

        if self.use_flanking_model:
            for f in ['FI', 'FD']:
                full_model.add_transition(full_model.start, states['flankleft'][f], 1.0)
                full_model.add_transition(states['flankright'][f], full_model.end, 1.0)

        LOGGER.info(str(xlink) + " crosslink edges were made")
        LOGGER.info("Baking")
        full_model.bake(True, merge="None")
        LOGGER.info("Baked")
        # LOGGER.info("There are %d total states and %d total edges", len(full_model.s), len(full_model.edges))

        return full_model

    def get_target(self, name) -> NucleotideSequence:
        """Get the target object with the given name if it exists in our data.

        Raises key error otherwise."""

        return self.sources[self._source_name_index_dict[name]]

    def align_from_fastx(
            self, query_fastx, source_fastx
    ) -> ProfalignAlignmentResultContainer:
        """Align all source sequences to the query sequence in the given FASTX
        (FASTA / FASTQ) files.

        NOTE: Assumes `query_fastx` contains only one sequence.
              If `query_fastx` contains more than one sequence, only the first
              will be considered.
        """

        queries = ingest_fastx_file(query_fastx)
        sources = ingest_fastx_file(source_fastx)
        sourceKeys = None

        # For now we enforce using only one query (the first one in the file):
        query = queries[0]

        return self.align(query, sources, sourceKeys)

    def align(self, query, sources, sourceKeys) -> ProfalignAlignmentResultContainer:
        """Align all sequences in `targets` to the given query sequence."""

        # Set our query and source properties:
        self._query_and_source_dependent_setup(query, sources, sourceKeys)
        for i, source in enumerate(self.sources):
            self._source_name_index_dict[source.name] = i

        # Create our panel here:
        panel = {self.query.name: self.query.sequence}
        for source in self.sources:
            panel[source.name] = source.sequence

        LOGGER.info("starting viterbi")
        logp, path = self.model.viterbi(query.sequence)
        LOGGER.info("finished viterbi")

        self.logp = logp

        ppath = []
        for p, (idx, state) in enumerate(path[1:-1]):
            # if (
            #         "start" not in state.name
            #         and ":RD" not in state.name
            #         and ":D" not in state.name
            # ):
            # print(state.name)
            ppath.append(f'{re.split(":", state.name)[0]}')

        self.editTrack: str = ""
        self.path: List[TesseraeAlignmentResult] = []

        self.__render(path, query, sources)

        self.alignment_result = ProfalignAlignmentResultContainer(self, path, self.path, None, logp)

        return self.alignment_result

    def __render(self, ppath, query, sources) -> None:
        """Trace back paths in the graph to create results.

        # todo: refactor this to return results with return statement
        Results are stored in self.path
        """

        # Prepare target sequence
        seqs = {}
        for s in range(len(sources)):
            seqs[sources[s].name] = sources[s].sequence
        sb = []
        pos_start = -1
        pos_end = -1
        pos_target = 1
        for p, (idx, state) in enumerate(ppath[1:-1]):
            if "start" not in state.name and "end" not in state.name and "RC" not in state.name and "FD" not in state.name:
                source_name, state = re.split(":", state.name)

                if state.startswith("D"):
                    sb.append("-")
                else:
                    if pos_start == -1:
                        pos_start = pos_target - 1

                    pos_end = pos_target - 1

                    sb.append(query.sequence[pos_target - 1])
                    pos_target += 1

        self.path.append(
            TesseraeAlignmentResult(query.name, "".join(sb), pos_start, pos_end)
        )

        # Prepare matching track
        sb = []
        pos_target = 1
        for p, (idx, state) in enumerate(ppath[1:-1]):
            if "start" not in state.name and "end" not in state.name and "FD" not in state.name:
                source_name, state = re.split(":", state.name)

                if source_name == "flankleft" or source_name == "flankright":
                    sb.append("*")
                elif state.startswith("M"):
                    ## TODO make this check for mismatches
                    sb.append("|")
                elif state.startswith("I"):
                    sb.append("^")
                else:
                    sb.append("~")
        self.editTrack = "".join(sb)

        # Prepare copying tracks
        current_track = "flankleft"
        sb = []
        pos_start = -1
        pos_end = -1
        total_bases = 0
        uppercase = True
        for p, (idx, state) in enumerate(ppath[1:-1]):
            if "start" not in state.name and "end" not in state.name and "FD" not in state.name:
                source_name, state = re.split(":", state.name)

                ## we have jumped to another sequence
                if current_track != source_name:
                    if (current_track != "flankleft" or self.use_flanking_model):
                        self.path.append(
                            TesseraeAlignmentResult(
                                current_track, "".join(sb), pos_start, pos_end
                            )
                        )
                    uppercase = True
                    pos_start = -1
                    pos_end = -1

                    current_track = source_name
                    sb = [repeat(" ", total_bases)]

                if "RC" not in state:
                    if source_name == "flankleft" or source_name == "flankright":
                        pos_start = 0
                        pos_end = pos_end + 1
                        sb.append(".")
                        total_bases += 1

                    elif state.startswith("I"):
                        total_bases += 1
                        sb.append("-")
                    else:
                        if not state.startswith("F"):
                            if "*" in state:
                                num = int(state[2:]) - 1
                            else:
                                num = int(state[1:]) - 1

                            if pos_start == -1:
                                pos_start = num
                            pos_end = num

                            c = seqs[source_name][num]
                            c = c.upper() if uppercase else c.lower()
                        # the state started with F, thus its an insetion so its pos_start/end are not valid numbers
                        else:
                            c = " "

                        total_bases += 1
                        sb.append(c)
        self.path.append(
            TesseraeAlignmentResult(current_track, "".join(sb), pos_start, pos_end)
        )
