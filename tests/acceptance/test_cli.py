from __future__ import annotations
import pathlib

import pysam

import pytest
from tesserae.__main__ import tesserae_cli

from ..util import EXPECTED_TEST_RESULTS, TEST_RESOURCES_FOLDER


def assert_sam_files_are_equal(samfile1, samfile2):
    """Asserts that the two given sam files are equal by comparing the file contents.

    Will correctly perform the comparison on bam or sam files."""

    with pysam.AlignmentFile(samfile1, "r") as f1, pysam.AlignmentFile(
        samfile2, "r"
    ) as f2:
        f1_header = f1.header.as_dict()
        f2_header = f2.header.as_dict()

        if 'PG' in f1_header:
            del f1_header['PG']

        if 'PG' in f2_header:
            del f2_header['PG']

        assert f1_header == f2_header

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


def test_tesserae_alignment_with_cli(tmpdir, capfd):
    out_file_path = pathlib.Path(tmpdir) / "test_tesserae_alignment_with_cli.bam"

    args = [
        "tesserae",
        "align",
        str(TEST_RESOURCES_FOLDER / "target_sequences.fasta"),
        str(TEST_RESOURCES_FOLDER / "query_sequence.fasta"),
        "-o",
        str(out_file_path),
    ]

    with pytest.raises(SystemExit):
        tesserae_cli(args[1:])

    # Now that we've run our data, let's do some basic checks:
    assert out_file_path.exists()
    assert out_file_path.is_file()

    # Check the contents of our output file:
    assert_sam_files_are_equal(out_file_path, EXPECTED_TEST_RESULTS)
