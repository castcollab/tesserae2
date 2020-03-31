import pytest
from tesserae.sequence import Sequence


class TestSequence:
    def test_trivial_operands(self):

        bases = "TAGCATGC"
        name = "test sequence"

        s = Sequence(name, bases)

        assert s.name == name
        assert s.sequence == bases

        assert len(s) == len(bases)
        assert str(s) == "{" + f"{name}: {bases}" + "}"
        assert s == Sequence(name, bases)

    @pytest.mark.parametrize("bases", ["atgcX", "L:", "3", "AUGC~", ":gagagag"])
    def test_invalid_sequence_creation_fails(self, bases):

        with pytest.raises(ValueError) as val_error:
            _ = Sequence("test", bases)
