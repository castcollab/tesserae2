import pytest
from tesserae.nucleotide_sequence import NucleotideSequence


class TestSequence:
    def test_trivial_operands(self):

        bases = "TAGCATGC"
        name = "test sequence"

        s = NucleotideSequence(name, bases)

        assert s.name == name
        assert s.sequence == bases

        assert len(s) == len(bases)
        assert str(s) == "{" + f"{name}: {bases}" + "}"
        assert s == NucleotideSequence(name, bases)

    @pytest.mark.parametrize("bases", ["atgcX", "L:", "3", "AUGC~", ":gagagag"])
    def test_invalid_sequence_creation_fails(self, bases):

        with pytest.raises(ValueError):
            _ = NucleotideSequence("test", bases)
