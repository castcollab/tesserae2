import re
from pathlib import Path
import pysam

from tesserae import cli

TEST_RESOURCES_FOLDER = Path(__file__).absolute().parent.parent / "resources"

EXPECTED_TEST_RESULTS = TEST_RESOURCES_FOLDER / "integration_test_expected_output.bam"

EXPECTED_SAM = """
@HD VN:1.0
@SQ SN:THE_QUERY LN:19
THE_FIRST_TARGET_0_6   0 THE_QUERY 2 60 7M  * 0 7  GTAGGCG      !!!!!!!
THE_SECOND_TARGET_7_18 0 THE_QUERY 9 60 12M * 0 12 AGATGACGCCAT !!!!!!!!!!!!
""".lstrip()
EXPECTED_SAM = re.sub(r" +", "\t", EXPECTED_SAM)


class SamFileRepr:
    def __init__(self, header, reads):
        self.header = header
        self.reads = reads

    def __str__(self):
        return "\n".join(self.header + self.reads) + "\n"


def quals_to_ascii(quals):
    return [str(chr(q + 33)) for q in quals]


def read_to_str(read, samfile):
    """Convert a pysam read to a human-readable format"""
    fields = str(read).split("\t")[:10]
    if fields[6] == "-1":
        fields[6] = "*"
    if fields[7] == "-1":
        fields[7] = "0"
    fields[3] = str(int(fields[3]) + 1)
    fields[2] = samfile.get_reference_name(int(fields[2]))
    fields.append("".join(quals_to_ascii(read.query_qualities)))
    return "\t".join(fields)


def load_bam(bam_file):

    samfile = pysam.AlignmentFile(bam_file, "rb")
    header = str(samfile.header).rstrip().split("\n")
    reads = []
    for read in samfile.fetch(until_eof=True):
        reads.append(read_to_str(read, samfile))
    return SamFileRepr(header=header, reads=reads)


class TestCli:
    def test_align_creates_expected_bam_file(self, tmpdir, capsys):
        # given
        out_file = Path(tmpdir) / "test_tesserae_alignment_with_cli.bam"
        args = [
            "tesserae",
            "align",
            "-v",
            "-r",
            str(TEST_RESOURCES_FOLDER / "query_sequence.fasta"),
            "-s",
            str(TEST_RESOURCES_FOLDER / "target_sequences.fasta"),
            "-o",
            str(out_file),
        ]

        # when
        cli.main(args)

        # then
        actual = load_bam(out_file)
        assert out_file.exists()
        assert out_file.is_file()
        assert EXPECTED_SAM == str(actual)

    def test_align_creates_expected_sam_stdout(self, tmpdir, capfd):
        # given
        args = [
            "tesserae",
            "align",
            "-v",
            "-r",
            str(TEST_RESOURCES_FOLDER / "query_sequence.fasta"),
            "-s",
            str(TEST_RESOURCES_FOLDER / "target_sequences.fasta"),
        ]

        # when
        cli.main(args)

        # then
        out, err = capfd.readouterr()
        assert EXPECTED_SAM == out
