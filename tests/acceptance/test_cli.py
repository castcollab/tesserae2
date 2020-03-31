import os
import pathlib
import filecmp
import pysam

from tesserae import cli

TEST_RESOURCES_FOLDER = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "resources"
)

EXPECTED_TEST_RESULTS = os.path.join(
    TEST_RESOURCES_FOLDER, "integration_test_expected_output.bam"
)


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
        assert filecmp.cmp(temp_file_path, EXPECTED_TEST_RESULTS)

        # TODO: Figure out why capsys fixture isn't working.
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
        # with pysam.AlignmentFile(
        #     tmp_stdout_file_path, "r"
        # ) as stdout_sam, pysam.AlignmentFile(
        #     EXPECTED_TEST_RESULTS, "rb"
        # ) as expected_bam:
        #
        #     assert stdout_sam.text == expected_bam.text
        #     assert stdout_sam.header.as_dict() == expected_bam.header.as_dict()
        #     assert stdout_sam.count() == expected_bam.count()
        #
        #     for actual, expected in zip(
        #         stdout_sam.fetch(until_eof=True), expected_bam.fetch(until_eof=True)
        #     ):
        #         assert actual == expected
