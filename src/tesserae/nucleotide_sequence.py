VALID_BASES = ("A", "a", "T", "t", "U", "u", "G", "g", "C", "c", "N", "n")


def validate_nucleic_acid_sequence(sequence: str, name: str = "") -> None:
    """Assert that if the given sequence is a string it contains only bases
    we expect from a nucleic acid sequence."""
    if len(name) != 0:
        name = "(" + name + ")"
    for i, base in enumerate(sequence):
        if base not in VALID_BASES:
            raise ValueError(
                f"Given sequence {name} contains unexpected base at "
                f"position {i}: {base}"
            )


class NucleotideSequence:
    """A named sequence of nucleotides."""

    def __init__(self, name: str, sequence: str) -> None:
        self._name = name
        self._sequence = sequence

        # validate_nucleic_acid_sequence(self._sequence, self._name)

    def __eq__(self, other) -> bool:
        if isinstance(other, NucleotideSequence):
            return (self.name == other.name) and (self.sequence == other.sequence)
        return False

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return "{" + f"{self.name}: {self.sequence}" + "}"

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._sequence
