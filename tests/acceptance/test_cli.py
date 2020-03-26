import os
import pathlib
import filecmp
import sys

from tesserae import cli

TEST_RESOURCES_FOLDER = os.path.join("tests", "resources")

EXPECTED_TEST_RESULTS = os.path.join("integration_test_expected_output.bam")


class TestCli:
    def test_tesserae_alignment_with_cli(self, tmpdir):
        """Tests that the CLI functions as expected and produces the expected results."""

        temp_file_path = os.path.join(tmpdir, "test_tesserae_alignment_with_cli.bam")

        sys.argv = [
            "tesserae",
            "-v",
            "align",
            "-r",
            os.path.join(TEST_RESOURCES_FOLDER, "query_sequence.fasta"),
            "-s",
            os.path.join(TEST_RESOURCES_FOLDER, "target_sequences.fasta"),
            "-o",
            temp_file_path,
        ]
        cli.main()

        # Now that we've run our data, let's do some basic checks:
        out_file_path = pathlib.Path(temp_file_path)
        assert out_file_path.exists() is True
        assert out_file_path.is_file() is True

        # OK, we have an output file.
        # Now let's check its contents:
        assert filecmp.cmp(temp_file_path, EXPECTED_TEST_RESULTS) is True
