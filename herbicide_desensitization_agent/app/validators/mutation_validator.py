import re

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
MUTATION = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")


def validate_mutation(mutation: str, sequence: str, protected: set[int]) -> tuple[str, int, str]:
    match = MUTATION.fullmatch(mutation)
    if not match:
        raise ValueError(f"Invalid single-substitution notation: {mutation}")
    source, raw_position, destination = match.groups()
    position = int(raw_position)
    if not 1 <= position <= len(sequence):
        raise ValueError(f"Mutation position outside sequence: {position}")
    if sequence[position - 1] != source:
        raise ValueError(f"Wild-type mismatch at {position}: expected {sequence[position - 1]}, got {source}")
    if source == destination:
        raise ValueError("Source and destination amino acids must differ")
    if position in protected:
        raise ValueError(f"Protected residue cannot be mutated: {position}")
    return source, position, destination

