VALID_BASES = ("A", "a", "T", "t", "U", "u", "G", "g", "C", "c", "N", "n")


def validate_nucleic_acid_sequence(sequence, name=""):
    """Assert that if the given sequence is a string it contains only bases
    we expect from a nucleic acid sequence."""
    assert isinstance(sequence, str), "Given sequence is not a string!"
    if len(name) != 0:
        name = "(" + name + ")"
    for i, base in enumerate(sequence):
        if base not in VALID_BASES:
            raise ValueError(
                f"Given sequence {name} contains unexpected base at "
                f"position {i}: {base}"
            )


class Sequence:
    """A named sequence of nucleotides."""

    def __init__(self, name, sequence):
        self._name = name
        self._sequence = sequence

        validate_nucleic_acid_sequence(self._sequence, self._name)

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return (self.name == other.name) and (self.sequence == other.sequence)
        return False

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return "{" + f"{self.name}: {self.sequence}" + "}"

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._sequence
