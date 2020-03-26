import os
import pathlib
import pysam
import pytest

from ..util import TEST_RESOURCES_FOLDER
from ..util import EXPECTED_TEST_RESULTS

from tesserae import cli


def assert_sam_files_are_equal(samfile1, samfile2):
    """Asserts that the two given sam files are equal by comparing the file contents.

    Will correctly perform the comparison on bam or sam files."""

    with pysam.AlignmentFile(samfile1, "r") as f1, pysam.AlignmentFile(
        samfile2, "r"
    ) as f2:
        assert f1.text == f2.text
        assert f1.header.as_dict() == f2.header.as_dict()

        f1_iter = f1.fetch(until_eof=True)
        f2_iter = f2.fetch(until_eof=True)

        for actual, expected in zip(f1_iter, f2_iter):
            assert actual == expected

        # Make sure both iterators are empty:
        with pytest.raises(StopIteration):
            _ = next(f1_iter)
        with pytest.raises(StopIteration):
            _ = next(f2_iter)


class TestCli:
    def test_tesserae_alignment_with_cli(self, tmpdir, capsys):
        temp_file_path = os.path.join(tmpdir, "test_tesserae_alignment_with_cli.bam")

        args = [
            "tesserae",
            "align",
            "-v",
            "-r",
            os.path.join(TEST_RESOURCES_FOLDER, "query_sequence.fasta"),
            "-s",
            os.path.join(TEST_RESOURCES_FOLDER, "target_sequences.fasta"),
            "-o",
            temp_file_path,
        ]
        cli.main(args)

        # Now that we've run our data, let's do some basic checks:
        out_file_path = pathlib.Path(temp_file_path)
        assert out_file_path.exists()
        assert out_file_path.is_file()

        # Check the contents of our output file:
        assert_sam_files_are_equal(temp_file_path, EXPECTED_TEST_RESULTS)

        # # TODO: Figure out why capsys fixture isn't working.
        # # Check the contents of stdout:
        # out, err = capsys.readouterr()
        #
        # # Create a sam file on disk for the stdout contents:
        # tmp_stdout_file_path = os.path.join(
        #     tmpdir, "test_tesserae_alignment_with_cli_stdout.sam"
        # )
        # with open(tmp_stdout_file_path, "w") as f:
        #     f.write(out)
        #
        # # Verify contents of our expected out and our sam file from stdout:
        # assert_sam_files_are_equal(tmp_stdout_file_path, EXPECTED_TEST_RESULTS)
