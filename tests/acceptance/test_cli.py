import pathlib
import subprocess

import pysam

import pytest
from tesserae import cli

from ..util import EXPECTED_TEST_RESULTS, TEST_RESOURCES_FOLDER


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

        for i, (actual, expected) in enumerate(zip(f1_iter, f2_iter)):
            assert actual == expected, f"Entry {i} is not equal in the two SAM files!"

        # Make sure both iterators are empty
        # (files have the same number of reads):
        # TODO: More elegant logging to indicate differing number of records:
        with pytest.raises(StopIteration):
            _ = next(f1_iter)
        with pytest.raises(StopIteration):
            _ = next(f2_iter)


class TestCli:
    @pytest.mark.parametrize("use_shell", (True, False))
    def test_tesserae_alignment_with_cli(self, tmpdir, capfd, use_shell):
        out_file_path = pathlib.Path(tmpdir) / "test_tesserae_alignment_with_cli.bam"

        args = [
            "tesserae",
            "-v",
            "align",
            "-r",
            str(TEST_RESOURCES_FOLDER / "query_sequence.fasta"),
            "-s",
            str(TEST_RESOURCES_FOLDER / "target_sequences.fasta"),
            "-o",
            str(out_file_path),
        ]

        # when
        if use_shell:
            subprocess.run(args)
        else:
            cli.main(args)

        # Now that we've run our data, let's do some basic checks:
        assert out_file_path.exists()
        assert out_file_path.is_file()

        # Check the contents of our output file:
        assert_sam_files_are_equal(out_file_path, EXPECTED_TEST_RESULTS)

        # Check the contents of stdout:
        out, _ = capfd.readouterr()

        # Create a sam file on disk for the stdout contents:
        tmp_stdout_file_path = (
            pathlib.Path(tmpdir) / "test_tesserae_alignment_with_cli_stdout.sam"
        )

        with open(tmp_stdout_file_path, "w") as f:
            f.write(out)

        # Verify contents of our expected out and our sam file from stdout:
        assert_sam_files_are_equal(tmp_stdout_file_path, EXPECTED_TEST_RESULTS)
