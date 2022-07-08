# distutils: language=c++
import sys
import itertools

import numpy

cimport cython
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint16_t
from libc.math cimport abs as cabs

from tesserae.utils import RefInterval
from tesserae.cpp_algorithm cimport transform, max_element, distance as iter_dist, find

NUCL = ['A', 'C', 'G', 'T']
NUCL_IX = {n: i for i, n in enumerate(NUCL)}
NUCL_TRANSITIONS = {
    ('A', 'G'), ('G', 'A'),  # Purines
    ('C', 'T'), ('T', 'C')  # Pyrimidines
}

DEF START = 0
DEF JUMP = 1

cdef float NINF = float('-inf')

# Type for elements of the pointer matrix
pointer_dtype = numpy.dtype([
    ('state', numpy.uint16),
    ('j', numpy.uint16)
])

ctypedef uint16_t state_t
ctypedef uint16_t pos_t

cdef packed struct Pointer:
    state_t state
    pos_t j


cdef unsigned char nucl_to_ix(char nucl):
    if nucl == b'A' or nucl == b'a':
        return 0
    elif nucl == b'C' or nucl == b'c':
        return 1
    elif nucl == b'G' or nucl == b'g':
        return 2
    elif nucl == b'T' or nucl == b't':
        return 3
    else:
        return 0xFF


cdef size_t cpp_argmax(vector[cython.numeric] to_search, cython.numeric value):
    return cabs(iter_dist(to_search.begin(), find(to_search.begin(), to_search.end(), value)))


cdef class TesseraeFactory:
    cdef float _pi_m  # Prob of starting in a match state
    cdef float _pi_i  # Prob of starting in an insert state
    cdef public float tau  # Prob of termination

    cdef public float delta  # gap initiation
    cdef public float eps  # gap extension
    cdef public float rho  # Jump between references

    cdef float _e_match
    cdef float _e_ts
    cdef float _e_tv

    cdef public object emissions_indel

    def __init__(self):
        # Transition matrix parameters
        self._pi_m = 0.75
        self._pi_i = 0.25
        self.tau = 0.001
        self.delta = 0.025
        self.eps = 0.75
        self.rho = 0.0001

        # Emission distribution parameters
        self._e_match = 0.9
        self._e_ts = 0.05
        self._e_tv = 0.025

        self.emissions_indel = numpy.full((4,), 0.25, dtype=numpy.float32)

    cpdef configure_start_end(self, float pi_m, float tau):
        assert 0 < pi_m < 1.0
        assert 0 < tau < 1.0

        self._pi_m = pi_m
        self._pi_i = 1 - pi_m
        self.tau = tau

        return self

    cpdef configure_gaps(self, float delta, float eps):
        assert 0 < delta < 1.0
        assert 0 < eps < 1.0

        self.delta = delta
        self.eps = eps

        return self

    cpdef configure_ref_jump(self, float rho):
        assert 0 <= rho < 1.0
        self.rho = rho

        return self

    cpdef configure_emissions(self, float e_match, float e_ts, object e_indel=None):
        assert e_match + e_ts < 1.0
        self._e_match = e_match
        self._e_ts = e_ts
        self._e_tv = (1 - e_match - e_ts) / 2

        if e_indel is not None:
            self.emissions_indel = e_indel

        return self

    property pi_i:
        def __get__(self):
            return self._pi_i

        def __set__(self, v):
            self._pi_i = v
            self._pi_m = 1 - v

    property pi_m:
        def __get__(self):
            return self._pi_m

        def __set__(self, v):
            self._pi_m = v
            self._pi_i = 1 - v


    cpdef Tesserae build(self, refs):
        cdef Tesserae hmm = Tesserae(refs)

        hmm.build_transition_matrix(self)
        hmm.build_emission_matrix(self)

        return hmm


cdef class AlignmentResult:
    cdef readonly states
    cdef readonly object viterbi
    cdef readonly object pointers

    def __init__(self, states, viterbi, pointers):
        # Matrix axes: (query pos, reference pos, state)
        self.states = states
        self.viterbi = viterbi
        self.pointers = pointers

    def get_log_likelihood(self):
        return self.viterbi[-1, -1, -1]

    cpdef get_aligned_ref_intervals(self):
        cdef float[:, :, ::1] vview = self.viterbi
        cdef Pointer[:, :, ::1] pview = self.pointers

        cdef Pointer curr_pointer = pview[-1, -1, -1]
        cdef pos_t i = pview.shape[0] - 1

        pointers = []
        while curr_pointer.state != START:
            pointers.append(curr_pointer)

            curr_state = self.states[curr_pointer.state]
            curr_pointer = pview[i, curr_pointer.j, curr_pointer.state]

            if curr_state[-2:] in {":M", ":I"}:
                i -= 1

        pointers.reverse()

        ref_intervals = []
        curr_interval = RefInterval()
        curr_op = None
        op_count = 0
        i = 0
        for curr_pointer in pointers:
            ref_state = self.states[curr_pointer.state]
            ref_id = ref_state[:-2]
            op = ref_state[-1]

            # Continuation of current reference
            if ref_state[-2:] in {":M", ":D", ":I"}:
                if curr_interval.ref_id is None:
                    curr_interval.ref_id = ref_id
                    curr_interval.qry_start = i

                if curr_interval.ref_start is None and op == "M":
                    curr_interval.ref_start = curr_pointer.j

                if curr_op is None:
                    curr_op = op

                if curr_op == op:
                    # Same CIGAR op (match/insert/deletion) as previous
                    op_count += 1
                else:
                    curr_interval.cigar_ops.append((op_count, curr_op))
                    op_count = 1

                # Keep updating end position as we trace the alignment, add one because Python intervals are
                # half-open intervals
                curr_interval.ref_end = curr_pointer.j + 1
                curr_interval.qry_end = i + 1
                curr_op = op

                if op in {"M", "I"}:
                    i += 1
            else:
                # New reference, append current reference interval to list, and create a new one
                curr_interval.cigar_ops.append((op_count, curr_op))
                ref_intervals.append(curr_interval)

                curr_interval = RefInterval()
                op_count = 0
                curr_op = None

        # Add the last CIGAR op and interval
        curr_interval.cigar_ops.append((op_count, curr_op))
        ref_intervals.append(curr_interval)

        return ref_intervals


cdef class Tesserae:

    cdef readonly unsigned int num_refs
    cdef readonly vector[string] ref_sequences
    cdef readonly vector[vector[pos_t]] ref_jump_sites

    cdef readonly object states
    cdef readonly object state_to_ref
    cdef readonly state_t num_states
    cdef readonly state_t ref_state_offset
    cdef readonly state_t states_per_ref

    # Transition and emission matrices
    # All log probabilities
    cdef readonly object transitions
    cdef readonly object emissions_match
    cdef readonly object emissions_indel

    def __init__(self, refs):
        self.num_refs = len(refs)

        self.states = ["begin", "jump"]
        self.ref_state_offset = len(self.states)
        self.states_per_ref = 3

        max_supported_refs = int((65535 - self.ref_state_offset) / self.states_per_ref)
        if len(refs) > max_supported_refs:
            raise ValueError(f"Too many references given. Maximum number of supported references: {max_supported_refs}")

        for i, r in enumerate(refs):
            if len(r) > 65535:
                raise ValueError(f"Reference {r.metadata['id']} is too long! Length of references can be at "
                                 f"most 65535.")

            self.states.append(f"{r.metadata['id']}:M")
            self.states.append(f"{r.metadata['id']}:I")
            self.states.append(f"{r.metadata['id']}:D")

            # Copy reference information to C++ objects for fast access in the `align` method
            # Also transform nucleotides to an integer 0-3, for easy access to the right index in matrices
            self.ref_sequences.push_back(bytes(map(nucl_to_ix, str(r).encode('utf-8'))))

            if 'jump_sites' in r.positional_metadata:
                # Obtain positions where recombination can happen from metadata (and always include start position)
                r.positional_metadata['jump_sites'][0] = True
                recomb_sites = list(r.positional_metadata['jump_sites'].nonzero()[0])
                self.ref_jump_sites.push_back(recomb_sites)
            else:
                # All positions as potential recombination sites
                self.ref_jump_sites.push_back(list(range(0, len(r))))

        self.states.append("end")

        self.num_states = len(self.states)

        self.transitions = None
        self.emissions_match = None
        self.emissions_indel = None

    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef build_transition_matrix(self, TesseraeFactory factory):
        self.transitions = numpy.zeros((self.num_states, self.num_states), dtype=numpy.float32)

        cdef float[:, ::1] matrix_view = self.transitions
        cdef unsigned int total_recomb_sites = sum(s.size() for s in self.ref_jump_sites)
        cdef state_t m, match_m_ix, ins_m_ix, del_m_ix

        for m in range(self.num_refs):
            # Determine indices within matrix of reference m
            match_m_ix = self.ref_state_offset + (m * self.states_per_ref)
            ins_m_ix = match_m_ix + 1
            del_m_ix = match_m_ix + 2

            # Starting in a state
            matrix_view[START, match_m_ix] = factory._pi_m / total_recomb_sites
            matrix_view[START, ins_m_ix] = factory._pi_i / total_recomb_sites

            # Match state to other states
            matrix_view[match_m_ix, match_m_ix] = 1.0 - (2*factory.delta) - factory.rho - factory.tau
            matrix_view[match_m_ix, ins_m_ix] = factory.delta
            matrix_view[match_m_ix, del_m_ix] = factory.delta
            matrix_view[match_m_ix, JUMP] = factory.rho
            matrix_view[match_m_ix, -1] = factory.tau

            # Insert state to other states
            matrix_view[ins_m_ix, match_m_ix] = 1.0 - factory.eps - factory.rho - factory.tau
            matrix_view[ins_m_ix, ins_m_ix] = factory.eps
            matrix_view[ins_m_ix, JUMP] = factory.rho
            matrix_view[ins_m_ix, -1] = factory.tau

            # Deletion state
            matrix_view[del_m_ix, match_m_ix] = 1.0 - factory.eps
            matrix_view[del_m_ix, del_m_ix] = factory.eps

            # Jump state
            matrix_view[JUMP, match_m_ix] = factory._pi_m / total_recomb_sites
            matrix_view[JUMP, ins_m_ix] = factory._pi_i / total_recomb_sites

        ix = (self.transitions > 0.0) & (self.transitions < 1.0)
        self.transitions[ix] = numpy.log(self.transitions[ix])
        self.transitions[self.transitions == 1.0] = 0
        self.transitions[self.transitions == 0.0] = NINF

    cdef build_emission_matrix(self, TesseraeFactory factory):
        self.emissions_match = numpy.zeros((4, 4), dtype=numpy.float32)
        for n1, n2 in itertools.combinations_with_replacement(NUCL, 2):
            ix1 = NUCL_IX[n1]
            ix2 = NUCL_IX[n2]

            if n1 == n2:
                self.emissions_match[ix1, ix2] = factory._e_match
                self.emissions_match[ix2, ix1] = factory._e_match
            elif (n1, n2) in NUCL_TRANSITIONS:
                self.emissions_match[ix1, ix2] = factory._e_ts
                self.emissions_match[ix2, ix1] = factory._e_ts
            else:
                self.emissions_match[ix1, ix2] = factory._e_tv
                self.emissions_match[ix2, ix1] = factory._e_tv

        self.emissions_match = numpy.log(self.emissions_match)
        self.emissions_indel = numpy.log(factory.emissions_indel)

    cpdef align(self, string query):
        """
        Use Viterbi algorithm to find the most likely alignment of the query to the panel of references.

        Parameters
        ----------
        query

        Returns
        -------
        """

        # Transform query nucleotides to integers 0-3 for easy emission matrix access
        transform(query.begin(), query.end(), query.begin(), nucl_to_ix)

        if query.size() > 65535:
            raise ValueError("Query sequence is too long! Can be at most 65535 long.")

        cdef pos_t i, j
        for i in range(query.size()):
            if query[i] >= 4:
                raise ValueError(f"Query sequence contains characters other than ACGT at position {i}.")

        cdef pos_t query_length = query.size()
        cdef pos_t max_ref_length = 0
        cdef string seq

        for seq in self.ref_sequences:
            max_ref_length = max(max_ref_length, seq.size())

        viterbi = numpy.full((query_length, max_ref_length, self.num_states), float('-inf'), dtype=numpy.float32)
        pointers = numpy.zeros((query_length, max_ref_length, self.num_states), dtype=pointer_dtype)
        cdef float[:, :, ::1] vview = viterbi
        cdef Pointer[:, :, ::1] pview = pointers

        # Global start state
        vview[0, 0, START] = 0  # log p, thus p = 1

        # Create views to transmission matrix and emission distribution matrix for fast access
        cdef float[:, ::1] tmatrix = self.transitions
        cdef float[:, ::1] emission_match_view = self.emissions_match
        cdef float[::1] emission_indel_view = self.emissions_indel

        # Set up base cases for match and insert states for each reference (can't start in a deletion state)
        cdef state_t ref, ref_m_ix, ref_i_ix, ref_d_ix
        for ref in range(self.num_refs):
            # State indices for the corresponding reference
            ref_m_ix = self.ref_state_offset + (ref * self.states_per_ref)
            ref_i_ix = ref_m_ix + 1

            # The alignment is allowed to start at any potential jump site
            for j in self.ref_jump_sites[ref]:
                # Start -> Match current ref at pos j
                vview[0, j, ref_m_ix] = (emission_match_view[query[0], self.ref_sequences[ref][j]] +
                                         tmatrix[START, ref_m_ix])
                pview[0, j, ref_m_ix] = Pointer(START, 0)

                # Start -> Insert vs ref at pos j
                vview[0, j, ref_i_ix] = emission_indel_view[query[0]] + tmatrix[START, ref_i_ix]
                pview[0, j, ref_i_ix] = Pointer(START, 0)

        # Base case for the jump state
        cdef float max_v_m = NINF
        cdef float max_v_i = NINF
        cdef Pointer max_m_pointer = Pointer(START, 0)
        cdef Pointer max_i_pointer = Pointer(START, 0)
        for ref in range(self.num_refs):
            # State indices for the corresponding reference
            ref_m_ix = self.ref_state_offset + (ref * self.states_per_ref)
            ref_i_ix = ref_m_ix + 1

            for j in self.ref_jump_sites[ref]:
                if max_v_m < vview[0, j, ref_m_ix]:
                    max_v_m = vview[0, j, ref_m_ix]
                    max_m_pointer = Pointer(ref_m_ix, j)

                if max_v_i < vview[0, j, ref_i_ix]:
                    max_v_i = vview[0, j, ref_i_ix]
                    max_i_pointer = Pointer(ref_i_ix, j)

        # Add transition probabilities (remember, log prob), and no emission prob because it's a silent state
        max_v_m += tmatrix[max_m_pointer.state, JUMP]
        max_v_i += tmatrix[max_i_pointer.state, JUMP]

        if max_v_m >= max_v_i:
            vview[0, 0, JUMP] = max_v_m
            pview[0, 0, JUMP] = max_m_pointer
        else:
            vview[0, 0, JUMP] = max_v_i
            pview[0, 0, JUMP] = max_i_pointer

        # All base cases initialized!
        # Prepare for the full recursion. Initialize C++ vectors that will contain the previous v-matrix values at
        # each iteration. By initializing them here we keep the inner loop below tight and fast.
        cdef vector[Pointer] prev_states_m = vector[Pointer](4)
        cdef vector[float] prev_v_m = vector[float](4)

        cdef vector[Pointer] prev_states_i = vector[Pointer](3)
        cdef vector[float] prev_v_i = vector[float](3)

        cdef vector[Pointer] prev_states_d = vector[Pointer](2)
        cdef vector[float] prev_v_d = vector[float](2)

        cdef vector[float].iterator max_prev
        cdef state_t argmax

        # Used to calculate jump state v-matrix values
        cdef float max_prev_m, max_prev_i, max_prev_r
        cdef Pointer pointer, pointer_m, pointer_i

        for i in range(1, query.size()):
            # First v-matrix values for non-silent states (i.e., reference match, insert or delete states)
            for j in range(0, max_ref_length):
                for ref in range(self.num_refs):
                    if j >= self.ref_sequences[ref].size():
                        continue

                    ref_m_ix = self.ref_state_offset + (ref * self.states_per_ref)
                    ref_i_ix = ref_m_ix + 1
                    ref_d_ix = ref_m_ix + 2

                    # (i, j) is in a match state, list potential previous states
                    prev_states_m.clear()
                    if j > 0:
                        prev_states_m.push_back(Pointer(ref_m_ix, j-1))
                        prev_states_m.push_back(Pointer(ref_i_ix, j-1))
                        prev_states_m.push_back(Pointer(ref_d_ix, j-1))
                    else:
                        # Dummy, but we add them because otherwise Cython's vector.size() and the actual size doesn't
                        # properly work. The actual vmatrix value will be set to NINF below, and thus not considered
                        # in the max() evaluation.
                        prev_states_m.push_back(Pointer(START, 0))
                        prev_states_m.push_back(Pointer(START, 0))
                        prev_states_m.push_back(Pointer(START, 0))

                    # JUMP state is independent of `j`
                    prev_states_m.push_back(Pointer(JUMP, 0))

                    prev_v_m.reserve(prev_states_m.size())
                    for prev_ix in range(prev_states_m.size()):
                        pointer = prev_states_m[prev_ix]
                        if pointer.state != START:
                            prev_v_m[prev_ix] = (tmatrix[pointer.state, ref_m_ix] +
                                                  vview[i-1, pointer.j, pointer.state])
                        else:
                            prev_v_m[prev_ix] = NINF

                    max_prev = max_element(prev_v_m.begin(), prev_v_m.end())
                    argmax = iter_dist(prev_v_m.begin(), max_prev)

                    vview[i, j, ref_m_ix] = deref(max_prev) + emission_match_view[query[i], self.ref_sequences[ref][j]]
                    pview[i, j, ref_m_ix] = prev_states_m[argmax]

                    # (i, j) is in an insert state, list potential previous states
                    prev_states_i.clear()
                    prev_states_i.push_back(Pointer(ref_m_ix, j))
                    prev_states_i.push_back(Pointer(ref_i_ix, j))

                    # JUMP state is independent of `j`
                    prev_states_i.push_back(Pointer(JUMP, 0))

                    prev_v_i.reserve(prev_states_i.size())
                    for prev_ix in range(prev_states_i.size()):
                        pointer = prev_states_i[prev_ix]
                        prev_v_i[prev_ix] = (tmatrix[pointer.state, ref_i_ix] +
                                             vview[i-1, pointer.j, pointer.state])

                    max_prev = max_element(prev_v_i.begin(), prev_v_i.end())
                    argmax = iter_dist(prev_v_i.begin(), max_prev)

                    vview[i, j, ref_i_ix] = deref(max_prev) + emission_indel_view[query[i]]
                    pview[i, j, ref_i_ix] = prev_states_i[argmax]

                    # (i, j) is in a deletion state, list potential previous states
                    if j > 0:
                        prev_states_d.clear()
                        prev_states_d.push_back(Pointer(ref_m_ix, j-1))
                        prev_states_d.push_back(Pointer(ref_d_ix, j-1))
                        for prev_ix in range(prev_states_d.size()):
                            pointer = prev_states_d[prev_ix]
                            prev_v_d[prev_ix] = (tmatrix[pointer.state, ref_d_ix] +
                                                 vview[i, pointer.j, pointer.state])

                        max_prev = max_element(prev_v_d.begin(), prev_v_d.end())
                        argmax = iter_dist(prev_v_d.begin(), max_prev)

                        vview[i, j, ref_d_ix] = deref(max_prev) + emission_indel_view[self.ref_sequences[ref][j]]
                        # vview[i, j, ref_d_ix] = deref(max_prev)
                        pview[i, j, ref_d_ix] = prev_states_d[argmax]

            # Calculate v-matrix values for silent states (the jump state is the only silent state)
            max_prev_m = NINF
            max_prev_i = NINF
            for ref in range(self.num_refs):
                ref_m_ix = self.ref_state_offset + (ref * self.states_per_ref)
                ref_i_ix = ref_m_ix + 1

                # Transition to JUMP state can come from any potential jump site in any reference
                for j in self.ref_jump_sites[ref]:
                    if vview[i, j, ref_m_ix] >= max_prev_m:
                        max_prev_m = vview[i, j, ref_m_ix]
                        pointer_m = Pointer(ref_m_ix, j)

                    if vview[i, j, ref_i_ix] >= max_prev_i:
                        max_prev_i = vview[i, j, ref_i_ix]
                        pointer_i = Pointer(ref_i_ix, j)

            max_prev_m += tmatrix[pointer_m.state, JUMP]
            max_prev_i += tmatrix[pointer_i.state, JUMP]

            if max_prev_m >= max_prev_i:
                vview[i, 0, JUMP] = max_prev_m
                pview[i, 0, JUMP] = pointer_m
            else:
                vview[i, 0, JUMP] = max_prev_i
                pview[i, 0, JUMP] = pointer_i

        # Transitions to end state
        cdef float max_prev_v = NINF
        cdef Pointer term_pointer = Pointer(0, 0)
        for ref in range(self.num_refs):
            ref_m_ix = self.ref_state_offset + (ref * self.states_per_ref)
            ref_i_ix = ref_m_ix + 1

            # In principle could end anywhere on any ref
            for j in range(max_ref_length):
                if j >= self.ref_sequences[ref].size():
                    continue

                if vview[-1, j, ref_m_ix] > max_prev_v:
                    max_prev_v = vview[-1, j, ref_m_ix]
                    term_pointer = Pointer(ref_m_ix, j)

                if vview[-1, j, ref_i_ix] > max_prev_v:
                    max_prev_v = vview[-1, j, ref_i_ix]
                    term_pointer = Pointer(ref_i_ix, j)

        max_prev_v += tmatrix[term_pointer.state, -1]

        # Just fill the whole `j` row with the same value, because the end state doesn't care about the
        # (i, j) coordinates
        vview[-1, :, -1] = max_prev_v
        pview[-1, :, -1] = term_pointer

        return AlignmentResult(self.states, viterbi, pointers)
