from numpy import sqrt

import pytest
from tesserae import Tesserae
from tesserae.nucleotide_sequence import NucleotideSequence
from tesserae.profmodel import Tesserae2

from ..util import TEST_RESOURCES_FOLDER

TEST_QUERY = NucleotideSequence("THE_QUERY", "GTAGGCGAGATGACGCCAT")
TEST_TARGETS = [
    NucleotideSequence("THE_FIRST_TARGET", "GTAGGCGAGTCCCGTTTATA"),
    NucleotideSequence("THE_SECOND_TARGET", "CCACAGAAGATGACGCCATT"),
]


class TesseraeFixture:
    def __init__(self):
        self.query = TEST_QUERY
        self.targets = TEST_TARGETS
        # when
        self.tesserae = Tesserae(mem_limit=False)
        self.alignment_results = self.tesserae.align(self.query, self.targets)


@pytest.fixture
def tesserae_fixture():
    return TesseraeFixture()


def assert_alignment_contents_are_correct(alignment, tesserae_object):

    assert len(alignment) == 3

    assert alignment[0].alignment_string == "GTAGGCGAGATGACGCCAT"

    assert alignment[1].alignment_string == "GTAGGCG"
    assert alignment[1].target_start_index == 0
    assert alignment[1].target_end_index == 6

    assert alignment[2].alignment_string == "       AGATGACGCCAT"
    assert alignment[2].target_start_index == 7
    assert alignment[2].target_end_index == 18

    assert tesserae_object.llk >= -27
    assert tesserae_object.llk <= -26

def assert_Tesserae2_alignment_contents_are_correct(alignment, tesserae_object):

    assert len(alignment) == 5

    assert alignment[0].alignment_string == "GTAGGCGAGATGACGCCAT"

    assert alignment[2].alignment_string == "GTAGGCG"
    assert alignment[2].target_start_index == 0
    assert alignment[2].target_end_index == 6

    assert alignment[4].alignment_string == "       AGATGACGCCAT"
    assert alignment[4].target_start_index == 7
    assert alignment[4].target_end_index == 18

    assert tesserae_object.logp >= -26
    assert tesserae_object.logp <= -25


class TestTesserae:
    @pytest.mark.parametrize(
        "query, targets, expected_query_name, expected_target_names",
        [
            pytest.param(
                TEST_QUERY,
                TEST_TARGETS,
                TEST_QUERY.name,
                [t.name for t in TEST_TARGETS],
                id="query Sequence target Sequence",
            ),
            pytest.param(
                NucleotideSequence("query", TEST_QUERY.sequence),
                TEST_TARGETS,
                "query",
                [t.name for t in TEST_TARGETS],
                id="query default target Sequence",
            ),
            pytest.param(
                TEST_QUERY,
                [
                    NucleotideSequence("sequence0", TEST_TARGETS[0].sequence),
                    NucleotideSequence("sequence1", TEST_TARGETS[1].sequence),
                ],
                TEST_QUERY.name,
                ["sequence0", "sequence1"],
                id="query Sequence target default",
            ),
        ],
    )
    def test_mosaic_alignment_on_short_query(
        self, query, targets, expected_query_name, expected_target_names
    ):
        """Test a mosaic alignment on a short query with the different inputs expected
        from Tesserae.align."""
        # Create our alignment:
        t = Tesserae()
        alignment_results = t.align(query, targets)

        # Assert that the names are correct:
        assert alignment_results[0].seq_name == expected_query_name
        for i in range(1, len(alignment_results)):
            assert alignment_results[i].seq_name == expected_target_names[i - 1]

        assert_alignment_contents_are_correct(alignment_results, t)

    def test_tesserae_properties(self):
        """Test that our properties function correctly."""

        # Create our alignment:
        t = Tesserae()
        _ = t.align(TEST_QUERY, TEST_TARGETS)

        # Individual property checks:
        assert t.query == TEST_QUERY
        assert t.targets == TEST_TARGETS

    def test_mosaic_alignment_on_short_query_from_fasta_files(self):
        # Create our alignment:

        t = Tesserae()
        a = t.align_from_fastx(
            TEST_RESOURCES_FOLDER / "query_sequence.fasta",
            TEST_RESOURCES_FOLDER / "target_sequences.fasta",
        )

        # Assert that the names are correct:
        assert a[0].seq_name == TEST_QUERY.name
        assert a[1].seq_name == TEST_TARGETS[0].name
        assert a[2].seq_name == TEST_TARGETS[1].name

        assert_alignment_contents_are_correct(a, t)

    def test_mosaic_alignment_on_short_query_and_two_templates_reduced_memory_usage(
        self, tesserae_fixture
    ):
        """Test a mosaic alignment on a short query with some memory limitations to
        ensure the results are correct."""
        # when
        t = Tesserae(mem_limit=True)
        alignment_results = t.align(tesserae_fixture.query, tesserae_fixture.targets)
        max_mem_limit = sqrt(len(tesserae_fixture.query)) + 1

        # then
        assert alignment_results == tesserae_fixture.alignment_results

        assert t.llk >= -27
        assert t.llk <= -26

        assert t.traceback_limit <= max_mem_limit
        assert t.num_states_to_save <= max_mem_limit

class TestProfalign:
    @pytest.mark.parametrize(
        "query, targets, expected_query_name, expected_target_names",
        [
            pytest.param(
                TEST_QUERY,
                TEST_TARGETS,
                TEST_QUERY.name,
                [t.name for t in TEST_TARGETS],
                id="query Sequence target Sequence",
            ),
            pytest.param(
                NucleotideSequence("query", TEST_QUERY.sequence),
                TEST_TARGETS,
                "query",
                [t.name for t in TEST_TARGETS],
                id="query default target Sequence",
            ),
            pytest.param(
                TEST_QUERY,
                [
                    NucleotideSequence("sequence0", TEST_TARGETS[0].sequence),
                    NucleotideSequence("sequence1", TEST_TARGETS[1].sequence),
                ],
                TEST_QUERY.name,
                ["sequence0", "sequence1"],
                id="query Sequence target default",
            ),
        ],
    )
    def test_mosaic_alignment_on_short_query(
        self, query, targets, expected_query_name, expected_target_names
    ):
        """Test a mosaic alignment on a short query with the different inputs expected
        from Tesserae.align."""
        # Create our alignment:
        t = Tesserae2()
        alignment_results = t.align(query, targets, None)

        # Assert that the names are correct:
        assert alignment_results[0].seq_name == expected_query_name
        assert alignment_results[1].seq_name == "flankleft"
        for i in range(1, len(expected_target_names)):
            assert alignment_results[2*i].seq_name == expected_target_names[i - 1]

        assert_Tesserae2_alignment_contents_are_correct(alignment_results, t)

    def test_tesserae_properties(self):
        """Test that our properties function correctly."""

        # Create our alignment:
        t = Tesserae2()
        _ = t.align(TEST_QUERY, TEST_TARGETS, None)

        # Individual property checks:
        assert t.query == TEST_QUERY
        assert t.sources == TEST_TARGETS

    def test_mosaic_alignment_on_medium_query_from_fasta_files_subsampling_10(self):
        # Create our alignment:

        t = Tesserae2(subsample_model_recombination=10)
        a = t.align_from_fastx(
            TEST_RESOURCES_FOLDER / "queryNoFlank.fasta",
            TEST_RESOURCES_FOLDER / "sources.fasta",
        )

        # Assert that the names are correct:
        assert a[0].seq_name == "query"
        assert a[1].seq_name == "flankleft"
        assert a[2].seq_name == "source1"
        assert a[3].seq_name == "recombination"
        assert a[4].seq_name == "source2"

        assert len(a) == 5

        assert a[0].alignment_string == "TCTTCGATATATAAGAAAAAAAATCGTCTAG-----TTCATCGTCTATGGGAGAGATTCCCTCAGCTAA"

        assert a[2].alignment_string == 'TCTTCGATATATAAGAAAAAAAATCGTCTA'
        assert a[2].target_start_index == 0
        assert a[2].target_end_index == 29

        assert a[4].alignment_string == '                              -ATAGCTTCATCGTCTATGGGAGAGATTCCCTCAGCTAA'
        assert a[4].target_start_index == 30
        assert a[4].target_end_index == 67

    def test_mosaic_alignment_on_medium_query_from_fasta_files(self):
        # Create our alignment:

        t = Tesserae2()
        a = t.align_from_fastx(
            TEST_RESOURCES_FOLDER / "queryNoFlank.fasta",
            TEST_RESOURCES_FOLDER / "sources.fasta",
        )

        # Assert that the names are correct:
        assert a[0].seq_name == "query"
        assert a[1].seq_name == "flankleft"
        assert a[2].seq_name == "source1"
        assert a[3].seq_name == "recombination"
        assert a[4].seq_name == "source2"

        assert len(a) == 5

        assert a[0].alignment_string == "TCTTCGATATATAAGAAAAAAAATCGTCTAGTTCATCGTCTATGGGAGAGATTCCCTCAGCTAA"

        assert a[2].alignment_string == "TCTTCGATATATAAGAAAAAAAATCGTCTAG"
        assert a[2].target_start_index == 0
        assert a[2].target_end_index == 30

        assert a[4].alignment_string == "                               TTCATCGTCTATGGGAGAGATTCCCTCAGCTAA"
        assert a[4].target_start_index == 35
        assert a[4].target_end_index == 67

    def test_mosaic_alignment_on_medium_query_from_fasta_files_with_flank(self):
        # Create our alignment:

        t = Tesserae2()
        a = t.align_from_fastx(
            TEST_RESOURCES_FOLDER / "query.fasta",
            TEST_RESOURCES_FOLDER / "sources.fasta",
        )

        # Assert that the names are correct:
        assert a[0].seq_name == "query"
        assert a[1].seq_name == "flankleft"
        assert a[2].seq_name == "source1"
        assert a[3].seq_name == "recombination"
        assert a[4].seq_name == "source2"
        assert a[5].seq_name == "flankright"

        assert len(a) == 6

        assert a[0].alignment_string == "TTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAATCTTCGATATATAAGAAAAAAAATCGTCTAGTTCATCGTCTATGGGAGAGATTCCCTCAGCTAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

        assert a[1].alignment_string == '..................................'
        assert a[1].target_start_index == 0
        assert a[1].target_end_index == 33

        assert a[2].alignment_string == "                                  TCTTCGATATATAAGAAAAAAAATCGTCTAG"
        assert a[2].target_start_index == 0
        assert a[2].target_end_index == 30

        assert a[4].alignment_string == "                                                                 TTCATCGTCTATGGGAGAGATTCCCTCAGCTAA"
        assert a[4].target_start_index == 35
        assert a[4].target_end_index == 67


    def test_mosaic_alignment_on_short_query_from_fasta_files(self):
        # Create our alignment:

        t = Tesserae2()
        a = t.align_from_fastx(
            TEST_RESOURCES_FOLDER / "query_sequence.fasta",
            TEST_RESOURCES_FOLDER / "target_sequences.fasta",
        )

        # Assert that the names are correct:
        assert a[0].seq_name == TEST_QUERY.name
        assert a[1].seq_name == "flankleft"
        assert a[2].seq_name == TEST_TARGETS[0].name
        assert a[3].seq_name == "recombination"
        assert a[4].seq_name == TEST_TARGETS[1].name

        assert_Tesserae2_alignment_contents_are_correct(a, t)
