from tesserae import Tesserae
from numpy import sqrt
import pytest
import os

# High-level test data storage to reference throughout the rest of the test cases:
TEST_DATA = ["GTAGGCGAGATGACGCCAT", ["GTAGGCGAGTCCCGTTTATA", "CCACAGAAGATGACGCCATT"]]

# Resources folder for test files:
TEST_RESOURCES_FOLDER = os.path.join("tests", "resources")


class TesseraeFixture:
    def __init__(self):
        self.query = TEST_DATA[0]
        self.targets = TEST_DATA[1]
        # when
        self.t = Tesserae(mem_limit=False)
        self.p = self.t.align(self.query, self.targets)


@pytest.fixture
def tesserae_fixture():
    return TesseraeFixture()


def get_default_target_name_list():
    """Returns a list of default target names for the targets in `test_data`."""
    return ["template" + str(i) for i in range(0, len(TEST_DATA[1]))]


def get_default_target_dictionary():
    """Returns a dictionary for the target test data:

     target_name -> target_sequence
     """
    d = {}
    for i in range(0, len(TEST_DATA[1])):
        d["template" + str(i)] = TEST_DATA[1][i]

    return d


class TestTesserae:
    """Class containing tests for the Tesserae Class."""

    @staticmethod
    def _assert_alignment_contents_are_correct(alignment, tesserae_object):
        """Asserts that the contents of the given alignment and tesserae object are correct for our short query set."""

        # Make sure we have the right number of alignments:
        assert len(alignment) == 3

        # Assert the contents are correct:
        assert alignment[0][1] == "GTAGGCGAGATGACGCCAT"

        assert alignment[1][1] == "GTAGGCG"
        assert alignment[1][2] == 0
        assert alignment[1][3] == 6

        assert alignment[2][1] == "       AGATGACGCCAT"
        assert alignment[2][2] == 7
        assert alignment[2][3] == 18

        assert tesserae_object.llk >= -27
        assert tesserae_object.llk <= -26

    @staticmethod
    @pytest.mark.parametrize(
        "query, targets, expected_query_name, expected_target_names",
        [
            pytest.param(
                TEST_DATA[0],
                TEST_DATA[1],
                "query",
                get_default_target_name_list(),
                id="query list target list",
            ),
            pytest.param(
                {"THE_QUERY": TEST_DATA[0]},
                {
                    "THE_FIRST_TARGET": TEST_DATA[1][0],
                    "THE_SECOND_TARGET": TEST_DATA[1][1],
                },
                "THE_QUERY",
                ["THE_FIRST_TARGET", "THE_SECOND_TARGET"],
                id="query dict target dict",
            ),
            pytest.param(
                TEST_DATA[0],
                {
                    "THE_FIRST_TARGET": TEST_DATA[1][0],
                    "THE_SECOND_TARGET": TEST_DATA[1][1],
                },
                "query",
                ["THE_FIRST_TARGET", "THE_SECOND_TARGET"],
                id="query list target dict",
            ),
            pytest.param(
                {"THE_QUERY": TEST_DATA[0]},
                TEST_DATA[1],
                "THE_QUERY",
                get_default_target_name_list(),
                id="query dict target list",
            ),
        ],
    )
    def test_mosaic_alignment_on_short_query(
        query, targets, expected_query_name, expected_target_names
    ):
        """Test a mosaic alignment on a short query with the different inputs expected from Tesserae.align."""
        # Create our alignment:
        t = Tesserae()
        a = t.align(query, targets)

        # Assert that the names are correct:
        assert a[0][0] == expected_query_name
        for i in range(1, len(a)):
            assert a[i][0] == expected_target_names[i - 1]

        TestTesserae._assert_alignment_contents_are_correct(a, t)

    def test_tesserae_properties(self):
        """Test that our properties function correctly."""

        # Create our alignment:
        t = Tesserae()
        a = t.align(TEST_DATA[0], TEST_DATA[1])

        # Individual property checks:
        assert t.query == {"query": TEST_DATA[0]}
        assert t.targets == get_default_target_dictionary()

    @staticmethod
    @pytest.mark.parametrize(
        "query, targets, expected_error_regex",
        [
            pytest.param("X", "ATGCATGC", r"0: X",),
            pytest.param("ATGCX", "ATGCATGC", r"4: X",),
            pytest.param("ATGCATGC,ATGC", "ATGCATGC", r"8: ,",),
        ],
    )
    def test_mosaic_alignment_on_short_query_raises_error_on_bad_sequence(
        query, targets, expected_error_regex
    ):
        # Create our alignment:

        with pytest.raises(ValueError, match=expected_error_regex):
            # Create our alignment:
            t = Tesserae()
            a = t.align(query, targets)

    def test_mosaic_alignment_on_short_query_from_fasta_files(self):
        # Create our alignment:

        t = Tesserae()
        a = t.align_from_fastx(
            os.path.join(TEST_RESOURCES_FOLDER, "query_sequence.fasta"),
            os.path.join(TEST_RESOURCES_FOLDER, "target_sequences.fasta"),
        )

        # Assert that the names are correct:
        assert a[0][0] == "THE_QUERY"
        assert a[1][0] == "THE_FIRST_TARGET"
        assert a[2][0] == "THE_SECOND_TARGET"

        TestTesserae._assert_alignment_contents_are_correct(a, t)

    def test_mosaic_alignment_on_short_query_and_two_templates_reduced_memory_usage(
        self, tesserae_fixture
    ):
        """Test a mosaic alignment on a short query with some memory limitations to ensure the results are correct."""
        # when
        t = Tesserae(mem_limit=True)
        p = t.align(tesserae_fixture.query, tesserae_fixture.targets)
        max_mem_limit = sqrt(len(tesserae_fixture.query)) + 1

        # then
        assert p == tesserae_fixture.p

        assert t.llk >= -27
        assert t.llk <= -26

        assert t.traceback_limit <= max_mem_limit
        assert t.states_to_save <= max_mem_limit
