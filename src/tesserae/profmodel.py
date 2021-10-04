import logging
import re

from dataclasses import dataclass
from typing import Dict, List, Sequence

import pysam
from pomegranate import *

from .nucleotide_sequence import NucleotideSequence

LOGGER = logging.getLogger(__name__)

# default transition probabilities
DEFAULT_DEL = 0.025
DEFAULT_EPS = 0.75
DEFAULT_REC = 0.0001
DEFAULT_TERM = 0.001

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
class TesseraeAlignmentResult:
    """
    Representation of the results of a Tesserae alignment
    """

    seq_name: str
    alignment_string: str
    target_start_index: int
    target_end_index: int


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
            self, threads=1, subsample_model_recombination=DEFAULT_SUBSAMPLING
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

        self._source_name_index_dict: Dict[str, int] = {}

    def _query_and_source_dependent_setup(self, query, sources, sourceKeys):
        """Perform all initialization that depends on the query or sources."""
        self.query = query
        self.sources = sources
        self.sourceKeys = sourceKeys
        self.model = self._make_full_model()

    def _make_global_alignment_model(self, source, keyArray):
        name = source.name
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
            extra = "*" if c+1 >= len(keyArray) or keyArray[c+1] else ""
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
        model.add_transition(model.start, s[f"{name}:D*1"], 0.05)
        model.add_transition(model.start, s[f"{name}:M*1"], 0.90)

        model.add_transition(s[f"{name}:I*0"], s[f"{name}:I*0"], 0.70)
        model.add_transition(s[f"{name}:I*0"], s[f"{name}:D*1"], 0.15)
        model.add_transition(s[f"{name}:I*0"], s[f"{name}:M*1"], 0.15)

        for c in range(1, len(source)):
            extra = "*" if keyArray[c] else ""
            extraPlus = "*" if c+1 >= len(keyArray) or keyArray[c+1] else ""
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
        full_model = self._make_flanking_model("flankleft")
        full_model.add_model(self._make_flanking_model("flankright"))

        for i in range(len(self.sources)):
            s = self.sources[i]
            if self.sourceKeys is None:
                source_key_to_use = self.subsample_model_recombination_sites(s, self.subsample_model_recombination)
            else:
                source_key_to_use = self.sourceKeys[i]

            assert len(source_key_to_use) == len(s), "Size mismatch for provided source recomination site array " \
                                                     "expected " +len(s) +" but found "+len(source_key_to_use)+" for " \
                                                                                                               "source:"+s.name

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
    ) -> List[TesseraeAlignmentResult]:
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

    def align(self, query, sources, sourceKeys) -> List[TesseraeAlignmentResult]:
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
        # print(path)

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

        return self.path

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