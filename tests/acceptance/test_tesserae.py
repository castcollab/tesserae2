import random
from random import choice

import pytest
from numpy import sqrt

from tesserae import Tesserae


def random_dna_string(length):
    return "".join([choice("ACGT") for _ in range(length)])


def repeat(string_to_expand, length):
    return (string_to_expand * (length // len(string_to_expand) + 1))[:length]


@pytest.fixture
def samples():
    random.seed(1)

    samples = {}
    for total_len in [100, 200, 500]:
        partial_len = total_len // 2
        query = random_dna_string(total_len)
        targets = [
            "".join(query[0:partial_len]) + random_dna_string(partial_len),
            random_dna_string(partial_len) + "".join(query[partial_len:total_len]),
        ]
        samples[total_len] = query, targets
    return samples


class TestTesseraeAcceptance:
    @classmethod
    def assertions(cls, p, targets, total_len):
        partial_len = total_len // 2

        assert len(p) == 3

        assert p[1][0] == "template0"
        assert p[1][1] == "".join(targets[0][0:partial_len])
        assert p[1][2] == 0
        assert p[1][3] == partial_len - 1

        assert p[2][0] == "template1"
        assert p[2][1] == repeat(" ", partial_len) + targets[1][partial_len:total_len]
        assert p[2][2] == partial_len
        assert p[2][3] == total_len - 1

    def test_mosaic_alignment_on_medium_queries_and_two_templates(self, samples):
        # given
        for total_len, (query, targets) in samples.items():
            # when
            t = Tesserae(mem_limit=False)
            p = t.align(query, targets)

            # then
            self.assertions(p, targets, total_len)

    def test_mosaic_alignment_on_medium_queries_and_two_templates_reduced_memory_usage(
        self, samples
    ):
        # given
        for total_len, (query, targets) in samples.items():
            # when
            t = Tesserae(mem_limit=True)
            p = t.align(query, targets)

            max_mem_limit = sqrt(total_len) + 2

            # then
            self.assertions(p, targets, total_len)
            assert t.traceback_limit <= max_mem_limit
            assert t.states_to_save <= max_mem_limit
