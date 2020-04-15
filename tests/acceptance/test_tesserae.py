import random
from random import choice

from numpy import sqrt

import pytest
from tesserae import Tesserae
from tesserae.sequence import Sequence


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
        query = Sequence("query", random_dna_string(total_len))
        targets = [
            Sequence(
                "target0",
                "".join(query.sequence[0:partial_len]) + random_dna_string(partial_len),
            ),
            Sequence(
                "target1",
                random_dna_string(partial_len)
                + "".join(query.sequence[partial_len:total_len]),
            ),
        ]
        samples[total_len] = query, targets
    return samples


class TestTesseraeAcceptance:
    @classmethod
    def __assertions__(cls, results, targets, total_len):
        partial_len = total_len // 2

        assert len(results) == 3

        assert results[1].seq_name == "target0"
        assert results[1].alignment_string == "".join(
            targets[0].sequence[0:partial_len]
        )
        assert results[1].target_start_index == 0
        assert results[1].target_end_index == partial_len - 1

        assert results[2].seq_name == "target1"
        assert (
            results[2].alignment_string
            == repeat(" ", partial_len) + targets[1].sequence[partial_len:total_len]
        )
        assert results[2].target_start_index == partial_len
        assert results[2].target_end_index == total_len - 1

    def test_mosaic_alignment_on_medium_queries_and_two_templates(self, samples):
        # given
        for total_len, (query, targets) in samples.items():
            # when
            t = Tesserae(mem_limit=False)
            alignment_results = t.align(query, targets)

            # then
            self.__assertions__(alignment_results, targets, total_len)

    def test_mosaic_alignment_on_medium_queries_and_two_templates_reduced_memory_usage(
        self, samples
    ):
        # given
        for total_len, (query, targets) in samples.items():
            # when
            t = Tesserae(mem_limit=True)
            target_alignment_results = t.align(query, targets)

            max_mem_limit = sqrt(total_len) + 2

            # then
            self.__assertions__(target_alignment_results, targets, total_len)
            assert t.traceback_limit <= max_mem_limit
            assert t.states_to_save <= max_mem_limit
