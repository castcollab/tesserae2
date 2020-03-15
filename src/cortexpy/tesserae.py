import numpy as np

# constants
SMALL = -1e32
STATES = 5

# default transition probabilities
DEFAULT_DEL = 0.025
DEFAULT_EPS = 0.75
DEFAULT_REC = 0.0001
DEFAULT_TERM = 0.001

# Nucleotide to index conversion
convert = {
    'G': 4, 'g': 4,
    'A': 3, 'a': 3,
    'C': 2, 'c': 2,
    'T': 1, 't': 1, 'U': 1, 'u': 1,
    'N': 0, 'n': 0,
}


def repeat(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]


def max_length(query, targets):
    maxl = len(query)

    for target in targets:
        maxl = max(maxl, len(target))

    return maxl


class Tesserae(object):
    pdel = DEFAULT_DEL
    peps = DEFAULT_EPS
    prho = DEFAULT_REC
    pterm = DEFAULT_TERM

    def __init__(
            self, mem_limit=False,
            pdel=DEFAULT_DEL, peps=DEFAULT_EPS,
            prho=DEFAULT_REC, pterm=DEFAULT_TERM
    ):
        self.pdel = pdel
        self.peps = peps
        self.prho = prho
        self.pterm = pterm

        self.ldel = np.log(self.pdel)
        self.leps = np.log(self.peps)
        self.lrho = np.log(self.prho)
        self.lterm = np.log(self.pterm)

        self.piM = 0.75
        self.piI = 1 - self.piM
        self.mm = 1 - 2 * self.pdel - self.prho - self.pterm
        self.gm = 1 - self.peps - self.prho - self.pterm
        self.dm = 1 - self.peps

        self.lpiM = np.log(self.piM)
        self.lpiI = np.log(self.piI)
        self.lmm = np.log(self.mm)
        self.lgm = np.log(self.gm)
        self.ldm = np.log(self.dm)

        self.emiss_gap_nt = [0.2, 0.2, 0.2, 0.2, 0.2]
        self.emiss_match_nt = [
            [0.2, 0.2, 0.2, 0.2, 0.2],
            [0.2, 0.9, 0.05, 0.025, 0.025],
            [0.2, 0.05, 0.9, 0.025, 0.025],
            [0.2, 0.025, 0.025, 0.9, 0.05],
            [0.2, 0.025, 0.025, 0.05, 0.9],
        ]

        self.editTrack = []
        self.mem_limit = mem_limit

    def __initialize(self, query, targets):
        self.nseq = len(targets)
        self.maxl = max_length(query, targets)
        self.qlen = len(query)

        if self.mem_limit:
            qlen_sqrt = np.sqrt(self.qlen)
            qsq_down = int(np.floor(qlen_sqrt))
            qsq_up = int(np.ceil(qlen_sqrt))
            # The traceback limit must be even!
            if qsq_up % 2 == 0:
                self.states_to_save = qsq_down
                self.traceback_limit = qsq_up
            else:
                self.states_to_save = qsq_up
                self.traceback_limit = qsq_down
            while self.states_to_save * self.traceback_limit < self.qlen:
                self.states_to_save += 1
        else:
            self.traceback_limit = self.qlen

        self.vt_m = np.full([2, self.nseq, self.maxl + 1], SMALL, dtype=np.float64)
        self.vt_i = np.full([2, self.nseq, self.maxl + 1], SMALL, dtype=np.float64)
        self.vt_d = np.full([2, self.nseq, self.maxl + 1], SMALL, dtype=np.float64)

        self.tb_m = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.maxl + 1],
            dtype=np.float64
        )
        self.tb_i = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.maxl + 1],
            dtype=np.float64
        )
        self.tb_d = np.zeros(
            [self.traceback_limit + 1, self.nseq, self.maxl + 1],
            dtype=np.float64
        )

        if self.mem_limit:
            self.saved_vt_m = np.full(
                [self.states_to_save + 1, self.nseq, self.maxl + 1],
                SMALL, dtype=np.float64
            )
            self.saved_vt_i = np.full(
                [self.states_to_save + 1, self.nseq, self.maxl + 1],
                SMALL, dtype=np.float64
            )
            self.saved_vt_d = np.full(
                [self.states_to_save + 1, self.nseq, self.maxl + 1],
                SMALL, dtype=np.float64
            )
            self.saved_states = []

        self.maxpath_copy = np.zeros([2 * self.maxl + 1], dtype=np.uint8)
        self.maxpath_state = np.zeros([2 * self.maxl + 1], dtype=np.uint8)
        self.maxpath_pos = np.zeros([2 * self.maxl + 1], dtype=np.int64)

        self.tmp = int(np.log(self.maxl) + 1)
        self.tb_divisor = np.power(10.0, self.tmp)

        self.combined_llk = 0.0

        self.sm = np.zeros([STATES, STATES], dtype=np.float64)
        self.lsm = np.zeros([STATES, STATES], dtype=np.float64)

        for i in range(0, STATES):
            for j in range(0, STATES):
                self.sm[i][j] = self.emiss_match_nt[i][j]
                self.lsm[i][j] = np.log(self.emiss_match_nt[i][j])

        self.si = np.zeros([STATES], dtype=np.float64)
        self.lsi = np.zeros([STATES], dtype=np.float64)

        for i in range(0, STATES):
            self.si[i] = self.emiss_gap_nt[i]
            self.lsi[i] = np.log(self.emiss_gap_nt[i])

        self.path = []

    def align(self, query, targets):
        self.__initialize(query, targets)

        panel = {"query": query}
        for i in range(0, len(targets)):
            panel[f'template{i}'] = targets[i]

        return self.__align_all(panel, targets)

    def __align_all(self, panel, targets):
        query = panel["query"]

        l1 = self.qlen
        size_l = 0.0
        for target in panel.values():
            size_l += len(target)

        size_l -= self.qlen
        lsize_l = np.log(size_l)

        max_r, pos_max, state_max, who_max = self.__initialization(query, targets, lsize_l)
        max_r, pos_max, state_max, who_max = self.__recurrence(
            query, targets, lsize_l, l1, max_r, pos_max,
            state_max, who_max, store_states=self.mem_limit
        )
        cp, pos_max, who_max, state_max = self.__termination(
            np.remainder(self.qlen, self.traceback_limit) if self.mem_limit else l1,
            pos_max, state_max, who_max, 2*self.maxl
        )
        l2 = self.traceback_limit
        if self.mem_limit and self.saved_states:
            self.saved_states.pop()
        while self.mem_limit and self.saved_states:
            max_r, pos_max_n, state_max_n, who_max_n, pos_target_n = self.saved_states.pop()
            idx = len(self.saved_states)
            self.vt_m[1 - (idx == 0)] = np.copy(self.saved_vt_m[idx])
            self.vt_i[1 - (idx == 0)] = np.copy(self.saved_vt_i[idx])
            self.vt_d[1 - (idx == 0)] = np.copy(self.saved_vt_d[idx])
            self.__recurrence(
                query, targets, lsize_l, l2, max_r, pos_max_n, state_max_n,
                who_max_n, offset=pos_target_n, l0=1 + (idx == 0), store_states=False
            )
            cp, pos_max, who_max, state_max = self.__termination(
                l2, pos_max, state_max, who_max, cp
            )
        self.__render(cp + 1, panel)

        return self.path

    def __render(self, cp, panel):
        # Prepare target sequence
        seqs = []
        for seqName in panel.keys():
            seqs.append((seqName, panel[seqName]))
        sb = []
        pos_start = -1
        pos_end = -1
        pos_target = 1
        for i in range(cp, 2 * self.maxl + 1):
            if self.maxpath_state[i] == 3:
                sb.append("-")
            else:
                if pos_start == -1:
                    pos_start = pos_target - 1

                pos_end = pos_target - 1

                sb.append(seqs[0][1][pos_target - 1])
                pos_target += 1
        self.path.append((seqs[0][0], "".join(sb), pos_start, pos_end))
        # Prepare matching track
        sb = []
        pos_target = 1
        for i in range(cp, 2 * self.maxl + 1):
            if self.maxpath_state[i] == 1:
                if seqs[0][1][pos_target - 1] == \
                   seqs[self.maxpath_copy[i] + 1][1][self.maxpath_pos[i] - 1]:
                    sb.append("|")
                else:
                    sb.append(" ")

                pos_target += 1
            elif self.maxpath_state[i] == 2:
                pos_target += 1
                sb.append("^")
            else:
                sb.append("~")
        self.editTrack = "".join(sb)
        # Prepare copying tracks
        current_track = seqs[self.maxpath_copy[cp] + 1][0]
        sb = []
        pos_start = -1
        pos_end = -1
        last_known_pos = -1
        uppercase = True
        for i in range(cp, 2 * self.maxl + 1):
            if i > cp and self.maxpath_copy[i] == self.maxpath_copy[i - 1] and \
                    np.abs(self.maxpath_pos[i] - self.maxpath_pos[i - 1]) > 1 or \
                    self.maxpath_pos[i] == last_known_pos + 1:
                self.path.append((current_track, "".join(sb), pos_start, pos_end))
                uppercase = not uppercase
                last_known_pos = self.maxpath_pos[i - 1]

                if pos_start != pos_end:
                    pos_start = self.maxpath_pos[i] - 1
                    pos_end = self.maxpath_pos[i] - 1

                current_track = seqs[self.maxpath_copy[i] + 1][0]
                sb = [repeat(' ', i - cp)]

            if i > cp and self.maxpath_copy[i] != self.maxpath_copy[i - 1]:
                self.path.append((current_track, "".join(sb), pos_start, pos_end))
                uppercase = True

                if pos_start != pos_end:
                    pos_start = self.maxpath_pos[i] - 1
                    pos_end = self.maxpath_pos[i] - 1

                current_track = seqs[self.maxpath_copy[i] + 1][0]
                sb = [repeat(' ', i - cp)]

            if self.maxpath_state[i] == 2:
                sb.append("-")
            else:
                c = seqs[self.maxpath_copy[i] + 1][1][self.maxpath_pos[i] - 1]
                c = c.upper() if uppercase else c.lower()

                if pos_start == -1:
                    pos_start = self.maxpath_pos[i] - 1

                pos_end = self.maxpath_pos[i] - 1

                sb.append(c)
        self.path.append((current_track, "".join(sb), pos_start, pos_end))

    def __to_traceback_indices(self, index):
        who = int(index / 10)
        state_float = index - who * 10
        state = int(state_float)
        pos = int((state_float - state) * self.tb_divisor + 1e-6)
        return who, state, pos

    def __termination(self, l1, pos_max, state_max, who_max, cp):
        self.maxpath_copy[cp] = who_max
        self.maxpath_state[cp] = state_max
        self.maxpath_pos[cp] = pos_max

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

            cp -= 1

            self.maxpath_copy[cp] = who_next
            self.maxpath_state[cp] = state_next
            self.maxpath_pos[cp] = pos_next

            who_max = who_next
            state_max = state_next
            pos_max = pos_next

            if self.maxpath_state[cp + 1] != 3:
                pos_target -= 1
        return cp, pos_max, who_max, state_max

    def __recurrence(
            self, query, targets, lsize_l, l1, max_r, pos_max,
            state_max, who_max, offset=0, l0=2, store_states=False
    ):
        who_max_n = 0
        state_max_n = 0
        pos_max_n = 0

        for pos_target in range(l0, l1 + 1):
            max_rn = SMALL + max_r
            seq = 0
            seq_10 = 0
            pos_target_trace = pos_target
            if self.mem_limit and store_states:
                pos_target_trace = pos_target % self.traceback_limit

            for target in targets:
                vt_m_base = max_r + self.lrho + self.lpiM - lsize_l
                vt_i_base = max_r + self.lrho + self.lpiI - lsize_l
                tb_base = who_max * 10 + state_max + pos_max / self.tb_divisor
                for pos_seq in range(1, len(target) + 1):
                    # Match
                    self.vt_m[1 - pos_target % 2][seq][pos_seq] = vt_m_base
                    self.tb_m[pos_target_trace][seq][pos_seq] = tb_base

                    vt_m_n, tb_m_n = max([
                        (self.vt_m[pos_target % 2][seq][pos_seq - 1] + self.lmm, 1),
                        (self.vt_i[pos_target % 2][seq][pos_seq - 1] + self.lgm, 2),
                        (self.vt_d[pos_target % 2][seq][pos_seq - 1] + self.ldm, 3)
                    ], key=lambda x: x[0])

                    if vt_m_n > self.vt_m[1 - pos_target % 2][seq][pos_seq]:
                        self.vt_m[1 - pos_target % 2][seq][pos_seq] = vt_m_n
                        self.tb_m[pos_target_trace][seq][pos_seq] = \
                            seq_10 + tb_m_n + (pos_seq - 1) / self.tb_divisor

                    # Add in state match
                    self.vt_m[1 - pos_target % 2][seq][pos_seq] += \
                        self.lsm[
                            convert[query[pos_target + offset - 1]]
                        ][
                            convert[target[pos_seq - 1]]
                        ]

                    # Insert
                    self.vt_i[1 - pos_target % 2][seq][pos_seq] = vt_i_base
                    self.tb_i[pos_target_trace][seq][pos_seq] = tb_base

                    vt_i_n, tb_i_n = max([
                        (self.vt_m[pos_target % 2][seq][pos_seq] + self.ldel, 1),
                        (self.vt_i[pos_target % 2][seq][pos_seq] + self.leps, 2)
                    ], key=lambda x: x[0])

                    if vt_i_n > self.vt_i[1 - pos_target % 2][seq][pos_seq]:
                        self.vt_i[1 - pos_target % 2][seq][pos_seq] = vt_i_n
                        self.tb_i[pos_target_trace][seq][pos_seq] = \
                            seq_10 + tb_i_n + pos_seq / self.tb_divisor

                    # Add in state insert
                    self.vt_i[1 - pos_target % 2][seq][pos_seq] += \
                        self.lsi[convert[query[pos_target + offset - 1]]]

                    # Delete
                    if (pos_target < l1 or (self.mem_limit and not store_states)) and pos_seq > 1:
                        vt_d_n, tb_d_n = max([
                            (self.vt_m[1 - pos_target % 2][seq][pos_seq - 1] + self.ldel, 1),
                            (self.vt_d[1 - pos_target % 2][seq][pos_seq - 1] + self.leps, 3)
                        ], key=lambda x: x[0])

                        self.vt_d[1 - pos_target % 2][seq][pos_seq] = vt_d_n
                        self.tb_d[pos_target_trace][seq][pos_seq] = \
                            seq_10 + tb_d_n + (pos_seq - 1) / self.tb_divisor

                    if self.vt_m[1 - pos_target % 2][seq][pos_seq] > max_rn:
                        max_rn = self.vt_m[1 - pos_target % 2][seq][pos_seq]
                        who_max_n = seq
                        state_max_n = 1
                        pos_max_n = pos_seq

                    if self.vt_i[1 - pos_target % 2][seq][pos_seq] > max_rn:
                        max_rn = self.vt_i[1 - pos_target % 2][seq][pos_seq]
                        who_max_n = seq
                        state_max_n = 2
                        pos_max_n = pos_seq

                seq += 1
                seq_10 += 10

            max_r = max_rn
            who_max = who_max_n
            state_max = state_max_n
            pos_max = pos_max_n
            if store_states and pos_target_trace == 0:
                idx = len(self.saved_states)
                self.saved_states.append((max_r, pos_max, state_max, who_max, pos_target))
                self.saved_vt_m[idx] = np.copy(self.vt_m[1])
                self.saved_vt_i[idx] = np.copy(self.vt_i[1])
                self.saved_vt_d[idx] = np.copy(self.vt_d[1])

        if (self.mem_limit and store_states) or not self.mem_limit:
            self.llk = max_r + self.lterm
            self.combined_llk += max_r + self.lterm

        return max_r, pos_max, state_max, who_max

    def __initialization(self, query, targets, lsize_l):
        who_max = 0
        state_max = 0
        pos_max = 0
        max_r = SMALL

        seq = 0
        for target in targets:
            for pos_seq in range(1, len(target) + 1):
                self.vt_m[0][seq][pos_seq] = self.lpiM - lsize_l + \
                    self.lsm[convert[query[0]]][convert[target[pos_seq - 1]]]
                self.vt_i[0][seq][pos_seq] = self.lpiI - lsize_l + \
                    self.lsi[convert[query[0]]]

                if pos_seq > 0:
                    vt_d_1, tb_d_1 = max([
                        (self.vt_m[0][seq][pos_seq - 1] + self.ldel, 1),
                        (self.vt_d[0][seq][pos_seq - 1] + self.leps, 3)
                    ], key=lambda x: x[0])
                    self.vt_d[0][seq][pos_seq] = vt_d_1
                    self.tb_d[0][seq][pos_seq] = seq * 10 + tb_d_1 + (pos_seq - 1) / self.tb_divisor

                if self.vt_m[0][seq][pos_seq] > max_r:
                    max_r = self.vt_m[0][seq][pos_seq]
                    who_max = seq
                    state_max = 1
                    pos_max = pos_seq

                if self.vt_i[0][seq][pos_seq] > max_r:
                    max_r = self.vt_i[0][seq][pos_seq]
                    who_max = seq
                    state_max = 2
                    pos_max = pos_seq

            seq += 1

        if self.mem_limit:
            self.saved_states.append((max_r, pos_max, state_max, who_max, 0))
            self.saved_vt_m[0] = np.copy(self.vt_m[0])
            self.saved_vt_i[0] = np.copy(self.vt_i[0])
            self.saved_vt_d[0] = np.copy(self.vt_d[0])

        return max_r, pos_max, state_max, who_max

    def __str__(self):
        max_name_length = 0

        for i in range(0, len(self.path)):
            name = "%s (%d-%d)" % (self.path[i][0], self.path[i][2], self.path[i][3])
            max_name_length = max(max_name_length, len(name))

        fmt = f'%{max_name_length}s'

        sb = ["\n"]

        for i in range(0, len(self.path)):
            name = "%s (%d-%d)" % (self.path[i][0], self.path[i][2], self.path[i][3])

            sb.append(fmt % name)
            sb.append(" ")
            sb.append(f'{self.path[i][1]}')
            sb.append("\n")

            if i == 0:
                sb.append(fmt % " ")
                sb.append(" ")
                sb.append(self.editTrack)
                sb.append("\n")

        sb.append("\n")
        sb.append(f'Mllk: {self.llk}')
        sb.append("\n")

        return "".join(sb)
