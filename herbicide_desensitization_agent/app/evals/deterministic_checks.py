from __future__ import annotations

from dataclasses import fields

from ..schemas.models import EvaluationPacket
from ..validators.mutation_validator import validate_mutation
from ..validators.provenance_validator import validate_provenance


def validate_evaluation_packet(
    packet: EvaluationPacket, sequence: str, protected_residues: set[int]
) -> list[str]:
    """Return deterministic validation failures; an empty list means pass."""
    failures: list[str] = []
    try:
        validate_mutation(packet.candidate.mutation, sequence, protected_residues)
    except ValueError as exc:
        failures.append(str(exc))
    try:
        validate_provenance(packet)
    except ValueError as exc:
        failures.append(str(exc))
    for score_field in fields(packet.scores):
        if not score_field.name.endswith("_score"):
            continue
        value = getattr(packet.scores, score_field.name)
        if not 0.0 <= value <= 1.0:
            failures.append(f"Score outside [0, 1]: {score_field.name}={value}")
    if packet.status == "APPROVED_FOR_ASSAY_PLANNING":
        failures.append("Prototype cannot autonomously approve candidates for assay planning")
    return failures
