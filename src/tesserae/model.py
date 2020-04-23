import collections
import itertools
import logging
import math
import sys
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pysam

from .nucleotide_sequence import NucleotideSequence

################################################################################

LOGGER = logging.getLogger(__name__)

################################################################################


# constants
SMALL = -1e32
STATES = 5

# default transition probabilities
DEFAULT_DEL = 0.025
DEFAULT_EPS = 0.75
DEFAULT_REC = 0.0001
DEFAULT_TERM = 0.001


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
class TesseraeAlignedSegment:
    """A segment of a Tesserae alignment."""

    seq_name: str
    alignment_string: str
    target_start_index: int
    target_end_index: int


@dataclass
class TesseraeAlignmentResult(collections.abc.Sequence):
    """Results of a Tesserae alignment."""

    path: List[TesseraeAlignedSegment]
    llk: float
    edit_track: str

    def __len__(self):
        return len(self.path)

    def __getitem__(self, item):
        return self.path[item]

    def __str__(self):
        max_name_length = 0

        for i in range(0, len(self.path)):
            name = "%s (%d-%d)" % (self.path[i][0], self.path[i][2], self.path[i][3],)
            max_name_length = max(max_name_length, len(name))

        fmt = f"%{max_name_length}s"

        sb = ["\n"]

        for i in range(0, len(self.path)):
            name = "%s (%d-%d)" % (self.path[i][0], self.path[i][2], self.path[i][3],)

            sb.append(fmt % name)
            sb.append(" ")
            sb.append(f"{self.path[i][1]}")
            sb.append("\n")

            if i == 0:
                sb.append(fmt % " ")
                sb.append(" ")
                sb.append(self.edit_track)
                sb.append("\n")

        sb.append("\n")
        sb.append(f"Mllk: {self.llk}")
        sb.append("\n")

        return "".join(sb)


def repeat(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]


def dump_query_and_targets_to_log(query, targets):
    """Dump the Sequence objects in query and targets to the log as a DEBUG message."""

    # Log our sequences:
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Query Sequence:")
        dump_sequence_to_log(query)
        LOGGER.debug("Target Sequences:")
        for target in targets:
            dump_sequence_to_log(target)


def dump_sequence_to_log(seq, spacing="    "):
    """Dump the given Sequence to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("%s%s -> %s", spacing, seq.name, seq.sequence)


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


@dataclass
class HMMIterationInfo:
    """These four values tend to be passed around between functions in the HMM a lot."""

    max_r: int
    pos_max: int
    state_max: int
    who_max: float


@dataclass
class RecurrenceTargetReturnVals:
    """Just something to avoid massive tuples"""

    vt_m_n: int
    vt_i_n: int
    vt_d_n: int
    tb_m_n: int
    tb_i_n: int
    tb_d_n: int


@dataclass
class TesseraeParameters:
    """A class for storing Tesserae model paramaters.

    An 'l' prefix indicates a log value.
    'lsm': emiss_match_nt
    'lsi': emiss_gap_nt
    """

    ldel: float
    leps: float
    lrho: float
    lterm: float
    lpiM: float
    lpiI: float
    lmm: float
    lgm: float
    ldm: float
    lsm: np.ndarray
    lsi: np.ndarray

    @classmethod
    def from_default_parms(cls):
        pdel = DEFAULT_DEL
        peps = DEFAULT_EPS
        prho = DEFAULT_REC
        pterm = DEFAULT_TERM
        piM = 0.75
        return cls.from_params(
            pdel=pdel,
            peps=peps,
            prho=prho,
            pterm=pterm,
            piM=piM,
            piI=1 - piM,
            mm=1 - 2 * pdel - prho - pterm,
            gm=1 - peps - prho - pterm,
            dm=1 - peps,
            emiss_gap_nt=np.array([0.2, 0.2, 0.2, 0.2, 0.2]),
            emiss_match_nt=np.array(
                [
                    [0.2, 0.2, 0.2, 0.2, 0.2],
                    [0.2, 0.9, 0.05, 0.025, 0.025],
                    [0.2, 0.05, 0.9, 0.025, 0.025],
                    [0.2, 0.025, 0.025, 0.9, 0.05],
                    [0.2, 0.025, 0.025, 0.05, 0.9],
                ]
            ),
        )

    @classmethod
    def from_params(
        cls, pdel, peps, prho, pterm, piM, piI, mm, gm, dm, emiss_gap_nt, emiss_match_nt
    ):
        """Builds a Tesserae model from regular-space (i.e. not log-space) parameters"""
        return cls(
            ldel=math.log(pdel),
            leps=math.log(peps),
            lrho=math.log(prho),
            lterm=math.log(pterm),
            lpiM=math.log(piM),
            lpiI=math.log(piI),
            lmm=math.log(mm),
            lgm=math.log(gm),
            ldm=math.log(dm),
            lsm=np.log(emiss_match_nt),
            lsi=np.log(emiss_gap_nt),
        )


class Tesserae:
    """Object to perform graph-based mosaic alignments between nucleic acid
    sequences.

    All parameters except for query and targets (at the moment) need to be
    specified at initialization time. In other words, all attributes should
    be considered read-only! Anything else may be madness.
    """

    def __init__(
        self, params: TesseraeParameters = None, mem_limit=False, threads=1,
    ):
        self.params = params or TesseraeParameters.from_default_parms()

        # Initialize the host of variables used elsewhere in Tesserae:
        self.num_states_to_save: int
        self.traceback_limit: int
        self.vt_m: np.ndarray[float]
        self.vt_i: np.ndarray[float]
        self.vt_d: np.ndarray[float]
        self.tb_m: np.ndarray[float]
        self.tb_i: np.ndarray[float]
        self.tb_d: np.ndarray[float]
        self.saved_vt_m: np.ndarray[float]
        self.saved_vt_i: np.ndarray[float]
        self.saved_vt_d: np.ndarray[float]
        self.saved_states: np.ndarray[float]
        self.maxpath_copy: np.ndarray[int]
        self.maxpath_state: np.ndarray[int]
        self.maxpath_pos: np.ndarray[int]
        self.tb_divisor: float
        self.llk = 0.0
        self.combined_llk = 0.0

        self.mem_limit: bool = mem_limit

        if sys.version_info < (3, 8) and threads > 1:
            raise NotImplementedError(
                "Threading is only supported with python3.8 and above!"
            )
        self.threads: int = threads

        # attributes specified in align()
        self.query: NucleotideSequence
        self.targets: Sequence[NucleotideSequence]
        self.ltarget_seq_size: float
        self.query_len: int
        self.max_seq_len: int
        self.nseq: int
        self._target_name_index_dict: Dict[str, int] = {}

    def _query_and_target_dependent_setup(self, query, targets):
        """Perform all initialization that depends on the query or targets."""
        self.query = query
        self.targets = targets
        self.ltarget_seq_size = math.log(sum(len(t.sequence) for t in self.targets))
        self.query_len = len(self.query)
        self.max_seq_len = max(
            len(s) for s in itertools.chain([self.query], self.targets)
        )
        self.nseq = len(self.targets)

        if self.mem_limit:
            qlen_sqrt = math.sqrt(self.query_len)
            qsq_down = math.floor(qlen_sqrt)
            qsq_up = math.ceil(qlen_sqrt)

            if qsq_up % 2 == 0:
                self.num_states_to_save = qsq_down
                self.traceback_limit = qsq_up
            else:
                self.num_states_to_save = qsq_up
                self.traceback_limit = qsq_down
            while self.num_states_to_save * self.traceback_limit < self.query_len:
                self.num_states_to_save += 1
        else:
            self.traceback_limit = self.query_len

        self.vt_m = np.full(
            [2, self.nseq, self.max_seq_len + 1], SMALL, dtype=np.float64
        )
        self.vt_i = np.full(
            [2, self.nseq, self.max_seq_len + 1], SMALL, dtype=np.float64
        )
        self.vt_d = np.full(
            [2, self.nseq, self.max_seq_len + 1], SMALL, dtype=np.float64
        )

        self.tb_m = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.max_seq_len + 1],
            dtype=np.float64,
        )
        self.tb_i = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.max_seq_len + 1],
            dtype=np.float64,
        )
        self.tb_d = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.max_seq_len + 1],
            dtype=np.float64,
        )

        if self.mem_limit:
            self.saved_vt_m = np.full(
                [self.num_states_to_save + 1, self.nseq, self.max_seq_len + 1],
                SMALL,
                dtype=np.float64,
            )
            self.saved_vt_i = np.full(
                [self.num_states_to_save + 1, self.nseq, self.max_seq_len + 1],
                SMALL,
                dtype=np.float64,
            )
            self.saved_vt_d = np.full(
                [self.num_states_to_save + 1, self.nseq, self.max_seq_len + 1],
                SMALL,
                dtype=np.float64,
            )
            self.saved_states = []

        self.maxpath_copy = np.zeros([2 * self.max_seq_len + 1], dtype=np.uint8)
        self.maxpath_state = np.zeros([2 * self.max_seq_len + 1], dtype=np.uint8)
        self.maxpath_pos = np.zeros([2 * self.max_seq_len + 1], dtype=np.int64)

        self.tb_divisor = np.power(10.0, int(math.log(self.max_seq_len) + 1))

    def get_target(self, name) -> NucleotideSequence:
        """Get the target object with the given name if it exists in our data.

        Raises key error otherwise."""

        return self.targets[self._target_name_index_dict[name]]

    def align_from_fastx(self, query_fastx, target_fastx) -> TesseraeAlignmentResult:
        """Align all target sequences to the query sequence in the given FASTX
        (FASTA / FASTQ) files.

        NOTE: Assumes `query_fastx` contains only one sequence.
              If `query_fastx` contains more than one sequence, only the first
              will be considered.
        """

        queries = ingest_fastx_file(query_fastx)
        targets = ingest_fastx_file(target_fastx)

        # For now we enforce using only one query (the first one in the file):
        query = queries[0]

        return self.align(query, targets)

    def align(self, query, targets) -> TesseraeAlignmentResult:
        """Align all sequences in `targets` to the given query sequence."""

        # Set our query and target properties:
        self._query_and_target_dependent_setup(query, targets)
        for i, target in enumerate(self.targets):
            self._target_name_index_dict[target.name] = i

        dump_query_and_targets_to_log(self.query, self.targets)

        iter_info = self.__initialization()
        iter_info = self.__recurrence(
            l1=self.query_len, iter_info=iter_info, store_states=self.mem_limit,
        )

        copy_position, pos_max, who_max, state_max = self.__termination(
            l1=self.query_len % self.traceback_limit
            if self.mem_limit
            else self.query_len,
            pos_max=iter_info.pos_max,
            state_max=iter_info.state_max,
            who_max=iter_info.who_max,
            copy_position=2 * self.max_seq_len,
        )

        if self.mem_limit and self.saved_states:
            self.saved_states.pop()

        while self.mem_limit and self.saved_states:
            iter_info, pos_target_n = self.saved_states.pop()

            idx = len(self.saved_states)
            self.vt_m[1 - (idx == 0)] = np.copy(self.saved_vt_m[idx])
            self.vt_i[1 - (idx == 0)] = np.copy(self.saved_vt_i[idx])
            self.vt_d[1 - (idx == 0)] = np.copy(self.saved_vt_d[idx])

            self.__recurrence(
                l1=self.traceback_limit,
                iter_info=iter_info,
                offset=pos_target_n,
                l0=1 + (idx == 0),
                store_states=False,
            )
            copy_position, pos_max, who_max, state_max = self.__termination(
                self.traceback_limit, pos_max, state_max, who_max, copy_position
            )

        return self.__render(copy_position + 1)

    def __render(self, copy_position: int) -> TesseraeAlignmentResult:
        """Trace back paths in the graph to create results."""
        builder = AlignmentTrackBuilder.from_tesserae(self, copy_position)
        path, edit_track = builder.build()
        return TesseraeAlignmentResult(path=path, llk=self.llk, edit_track=edit_track)

    def __to_traceback_indices(self, index):
        who = int(index / 10)
        state_float = index - who * 10
        state = int(state_float)
        pos = int((state_float - state) * self.tb_divisor + 1e-6)
        return who, state, pos

    def __termination(self, l1, pos_max, state_max, who_max, copy_position):
        self.maxpath_copy[copy_position] = who_max
        self.maxpath_state[copy_position] = state_max
        self.maxpath_pos[copy_position] = pos_max

        pos_target = l1
        who_next = 0
        state_next = 0
        pos_next = 0
        while pos_target >= 1:
            if state_max == 1:
                who_next, state_next, pos_next = self.__to_traceback_indices(
                    self.tb_m[pos_target][who_max][pos_max]
                )
            elif state_max == 2:
                who_next, state_next, pos_next = self.__to_traceback_indices(
                    self.tb_i[pos_target][who_max][pos_max]
                )
            elif state_max == 3:
                who_next, state_next, pos_next = self.__to_traceback_indices(
                    self.tb_d[pos_target][who_max][pos_max]
                )

            copy_position -= 1

            self.maxpath_copy[copy_position] = who_next
            self.maxpath_state[copy_position] = state_next
            self.maxpath_pos[copy_position] = pos_next

            who_max = who_next
            state_max = state_next
            pos_max = pos_next

            if self.maxpath_state[copy_position + 1] != 3:
                pos_target -= 1
        return copy_position, pos_max, who_max, state_max

    @staticmethod
    def recurrence_target(
        query,
        target,
        vt_m,
        vt_i,
        vt_d,
        vt_m_base,
        vt_i_base,
        tb_base,
        max_rn,
        seq,
        pos_target,
        offset,
        tb_divisor,
        params: TesseraeParameters,
    ):
        seq_10 = seq * 10

        targ_idx = CONVERT[query.sequence[pos_target + offset - 1]]

        em_m = np.insert(
            np.array(
                [
                    params.lsm[targ_idx][CONVERT[target.sequence[i]]]
                    if i < len(target)
                    else SMALL
                    for i in range(0, len(vt_m) - 1)
                ],
                dtype=np.float64,
            ),
            0,
            SMALL,
        )
        em_m_tb = np.insert(
            (seq_10 + np.arange(0, len(vt_m) - 1, dtype=np.float64) / tb_divisor),
            0,
            tb_base,
        )
        vt_m_mat = (
            np.c_[
                np.full(len(vt_m), vt_m_base, dtype=np.float64),
                np.roll(vt_m, 1) + params.lmm,
                np.roll(vt_i, 1) + params.lgm,
                np.roll(vt_d, 1) + params.ldm,
            ]
            + np.c_[em_m, em_m, em_m, em_m]
        )
        tb_m_mat = np.c_[
            np.full(len(vt_m), tb_base, dtype=np.float64),
            em_m_tb + 1,
            em_m_tb + 2,
            em_m_tb + 3,
        ]

        tb_m_n = vt_m_mat.argmax(1)
        index_selector = tb_m_n + np.arange(len(tb_m_n)) * 4
        vt_m_n = vt_m_mat.ravel()[index_selector]
        tb_m_n = tb_m_mat.ravel()[index_selector]

        em_i = np.insert(
            np.full(len(vt_i) - 1, params.lsi[targ_idx], dtype=np.float64), 0, SMALL
        )
        em_i_tb = seq_10 + np.arange(len(vt_m), dtype=np.float64) / tb_divisor
        vt_i_mat = (
            np.c_[
                np.full(len(vt_i), vt_i_base, dtype=np.float64),
                vt_m + params.ldel,
                vt_i + params.leps,
            ]
            + np.c_[em_i, em_i, em_i]
        )
        tb_i_mat = np.c_[
            np.full(len(vt_i), tb_base, dtype=np.float64), em_i_tb + 1, em_i_tb + 2,
        ]

        tb_i_n = vt_i_mat.argmax(1)
        index_selector = tb_i_n + np.arange(len(tb_i_n)) * 3
        vt_i_n = vt_i_mat.ravel()[index_selector]
        tb_i_n = tb_i_mat.ravel()[index_selector]

        idx_m = vt_m_n.argmax()
        idx_i = vt_i_n.argmax()

        if vt_m_n[idx_m] > max_rn:
            max_rn = vt_m_n[idx_m]
            who_max_n = seq
            state_max_n = 1
            pos_max_n = idx_m

        if vt_i_n[idx_i] > max_rn:
            max_rn = vt_i_n[idx_i]
            who_max_n = seq
            state_max_n = 2
            pos_max_n = idx_i

        # Use traditional python-list instead of numpy for
        # the delete-vector in order to exploit the best of two worlds
        vt_d_next = list(np.roll(vt_m_n, 1) + params.ldel)
        tb_d_next = [1] * len(vt_m_n)
        for pos_seq in range(2, len(target) + 1):
            vt_d_n = vt_d_next[pos_seq - 1] + params.leps
            if vt_d_next[pos_seq] < vt_d_n:
                vt_d_next[pos_seq] = vt_d_n
                tb_d_next[pos_seq] = 3

        vt_d_n = np.array(vt_d_next)
        tb_d_n = em_m_tb + tb_d_next

        return (
            HMMIterationInfo(
                max_r=max_rn,
                pos_max=pos_max_n,
                state_max=state_max_n,
                who_max=who_max_n,
            ),
            RecurrenceTargetReturnVals(vt_m_n, vt_i_n, vt_d_n, tb_m_n, tb_i_n, tb_d_n),
        )

    def __recurrence(
        self, l1, iter_info: HMMIterationInfo, offset=0, l0=2, store_states=False,
    ) -> HMMIterationInfo:
        # Use no more threads than sequences, since all
        # parallelizations are on the individual sequences.
        recurrence_threads = min(self.threads, self.nseq)

        for pos_target in range(l0, l1 + 1):
            max_rn = SMALL + iter_info.max_r
            pos_target_trace = pos_target
            if self.mem_limit and store_states:
                pos_target_trace = pos_target % self.traceback_limit

            vt_m_base = (
                iter_info.max_r
                + self.params.lrho
                + self.params.lpiM
                - self.ltarget_seq_size
            )
            vt_i_base = (
                iter_info.max_r
                + self.params.lrho
                + self.params.lpiI
                - self.ltarget_seq_size
            )
            tb_base = (
                iter_info.who_max * 10
                + iter_info.state_max
                + iter_info.pos_max / self.tb_divisor
            )

            argvec = []
            for target_idx, target in enumerate(self.targets):
                argvec.append(
                    (
                        self.query,
                        target,
                        self.vt_m[pos_target % 2][target_idx],
                        self.vt_i[pos_target % 2][target_idx],
                        self.vt_d[pos_target % 2][target_idx],
                        vt_m_base,
                        vt_i_base,
                        tb_base,
                        max_rn,
                        target_idx,
                        pos_target,
                        offset,
                        self.tb_divisor,
                        self.params,
                    )
                )

            # todo: check if the list comprehension and map are a memory hog.
            #  If so convert to generators and imap
            if recurrence_threads > 1:
                from multiprocessing.pool import ThreadPool as Pool

                with Pool(recurrence_threads) as pool:
                    rvec = pool.map(
                        lambda args: Tesserae.recurrence_target(*args), argvec
                    )
            else:
                rvec = [Tesserae.recurrence_target(*args) for args in argvec]
            iter_info, _ = max(rvec, key=lambda x: x[0].max_r)
            for target_idx, (_, ret_vals) in enumerate(rvec):
                self.vt_m[1 - pos_target % 2][target_idx] = ret_vals.vt_m_n
                self.vt_i[1 - pos_target % 2][target_idx] = ret_vals.vt_i_n
                self.vt_d[1 - pos_target % 2][target_idx] = ret_vals.vt_d_n
                self.tb_m[pos_target_trace][target_idx] = ret_vals.tb_m_n
                self.tb_i[pos_target_trace][target_idx] = ret_vals.tb_i_n
                self.tb_d[pos_target_trace][target_idx] = ret_vals.tb_d_n

            if store_states and pos_target_trace == 0:
                idx = len(self.saved_states)
                self.saved_states.append((iter_info, pos_target))
                self.saved_vt_m[idx] = np.copy(self.vt_m[1])
                self.saved_vt_i[idx] = np.copy(self.vt_i[1])
                self.saved_vt_d[idx] = np.copy(self.vt_d[1])

        if (self.mem_limit and store_states) or not self.mem_limit:
            self.llk = iter_info.max_r + self.params.lterm
            self.combined_llk += iter_info.max_r + self.params.lterm

        return iter_info

    def __initialization(self):
        who_max = 0
        state_max = 0
        pos_max = 0
        max_r = SMALL

        for target_idx, target in enumerate(self.targets):
            for target_pos in range(1, len(target) + 1):
                self.vt_m[0][target_idx][target_pos] = (
                    self.params.lpiM
                    - self.ltarget_seq_size
                    + self.params.lsm[CONVERT[self.query.sequence[0]]][
                        CONVERT[target.sequence[target_pos - 1]]
                    ]
                )
                self.vt_i[0][target_idx][target_pos] = (
                    self.params.lpiI
                    - self.ltarget_seq_size
                    + self.params.lsi[CONVERT[self.query.sequence[0]]]
                )

                if target_pos > 0:
                    vt_d_1, tb_d_1 = max(
                        [
                            (
                                self.vt_m[0][target_idx][target_pos - 1]
                                + self.params.ldel,
                                1,
                            ),
                            (
                                self.vt_d[0][target_idx][target_pos - 1]
                                + self.params.leps,
                                3,
                            ),
                        ],
                        key=lambda x: x[0],
                    )
                    self.vt_d[0][target_idx][target_pos] = vt_d_1
                    self.tb_d[0][target_idx][target_pos] = (
                        target_idx * 10 + tb_d_1 + (target_pos - 1) / self.tb_divisor
                    )

                if self.vt_m[0][target_idx][target_pos] > max_r:
                    max_r = self.vt_m[0][target_idx][target_pos]
                    who_max = target_idx
                    state_max = 1
                    pos_max = target_pos

                if self.vt_i[0][target_idx][target_pos] > max_r:
                    max_r = self.vt_i[0][target_idx][target_pos]
                    who_max = target_idx
                    state_max = 2
                    pos_max = target_pos

        if self.mem_limit:
            self.saved_states.append(
                (HMMIterationInfo(max_r, pos_max, state_max, who_max), 0)
            )
            self.saved_vt_m[0] = np.copy(self.vt_m[0])
            self.saved_vt_i[0] = np.copy(self.vt_i[0])
            self.saved_vt_d[0] = np.copy(self.vt_d[0])

        return HMMIterationInfo(max_r, pos_max, state_max, who_max)


@dataclass
class AlignmentTrackBuilder:
    seqs: List[Tuple[str, str]]
    max_seq_len: int
    maxpath_state: List[int]
    maxpath_copy: List[int]
    maxpath_pos: List[int]
    copy_position: int

    @classmethod
    def from_tesserae(cls, tr: Tesserae, copy_position):
        seqs = [(s.name, s.sequence) for s in itertools.chain([tr.query], tr.targets)]
        return cls(
            seqs=seqs,
            max_seq_len=tr.max_seq_len,
            maxpath_state=tr.maxpath_state,
            maxpath_copy=tr.maxpath_copy,
            maxpath_pos=tr.maxpath_pos,
            copy_position=copy_position,
        )

    def build(self):
        first_track = self._build_first_track()

        # Prepare matching track
        edit_track = self._build_edit_track()

        # Prepare copying tracks
        path = self._build_copying_tracks()
        return [first_track] + path, edit_track

    def _build_first_track(self) -> TesseraeAlignedSegment:
        string_builder = []
        pos_start = -1
        pos_end = -1
        pos_target = 1
        for i in range(self.copy_position, 2 * self.max_seq_len + 1):
            if self.maxpath_state[i] == 3:
                string_builder.append("-")
            else:
                if pos_start == -1:
                    pos_start = pos_target - 1

                pos_end = pos_target - 1

                string_builder.append(self.seqs[0][1][pos_target - 1])
                pos_target += 1
        return TesseraeAlignedSegment(
            self.seqs[0][0], "".join(string_builder), pos_start, pos_end
        )

    def _build_copying_tracks(self) -> List[TesseraeAlignedSegment]:
        current_track = self.seqs[self.maxpath_copy[self.copy_position] + 1][0]
        string_builder: List[str] = []
        pos_start = -1
        pos_end = -1
        last_known_pos = -1
        uppercase = True
        path = []
        for i in range(self.copy_position, 2 * self.max_seq_len + 1):
            if (
                i > self.copy_position
                and self.maxpath_copy[i] == self.maxpath_copy[i - 1]
                and np.abs(self.maxpath_pos[i] - self.maxpath_pos[i - 1]) > 1
                or self.maxpath_pos[i] == last_known_pos + 1
            ):
                path.append(
                    TesseraeAlignedSegment(
                        current_track, "".join(string_builder), pos_start, pos_end
                    )
                )
                uppercase = not uppercase
                last_known_pos = self.maxpath_pos[i - 1]

                if pos_start != pos_end:
                    pos_start = self.maxpath_pos[i] - 1
                    pos_end = self.maxpath_pos[i] - 1

                current_track = self.seqs[self.maxpath_copy[i] + 1][0]
                string_builder = [repeat(" ", i - self.copy_position)]

            if (
                i > self.copy_position
                and self.maxpath_copy[i] != self.maxpath_copy[i - 1]
            ):
                path.append(
                    TesseraeAlignedSegment(
                        current_track, "".join(string_builder), pos_start, pos_end
                    )
                )
                uppercase = True

                if pos_start != pos_end:
                    pos_start = self.maxpath_pos[i] - 1
                    pos_end = self.maxpath_pos[i] - 1

                current_track = self.seqs[self.maxpath_copy[i] + 1][0]
                string_builder = [repeat(" ", i - self.copy_position)]

            if self.maxpath_state[i] == 2:
                string_builder.append("-")
            else:
                c = self.seqs[self.maxpath_copy[i] + 1][1][self.maxpath_pos[i] - 1]
                c = c.upper() if uppercase else c.lower()

                if pos_start == -1:
                    pos_start = self.maxpath_pos[i] - 1

                pos_end = self.maxpath_pos[i] - 1

                string_builder.append(c)
        path.append(
            TesseraeAlignedSegment(
                current_track, "".join(string_builder), pos_start, pos_end
            )
        )
        return path

    def _build_edit_track(self) -> str:
        string_builder = []
        pos_target = 1
        for i in range(self.copy_position, 2 * self.max_seq_len + 1):
            if self.maxpath_state[i] == 1:
                if (
                    self.seqs[0][1][pos_target - 1]
                    == self.seqs[self.maxpath_copy[i] + 1][1][self.maxpath_pos[i] - 1]
                ):
                    string_builder.append("|")
                else:
                    string_builder.append(" ")

                pos_target += 1
            elif self.maxpath_state[i] == 2:
                pos_target += 1
                string_builder.append("^")
            else:
                string_builder.append("~")
        return "".join(string_builder)
