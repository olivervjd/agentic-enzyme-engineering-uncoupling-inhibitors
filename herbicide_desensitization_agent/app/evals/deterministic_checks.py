from __future__ import annotations

from dataclasses import fields

from ..schemas.models import EvaluationPacket, TargetRegistryEntry
from ..validators.mutation_validator import validate_mutation
from ..validators.provenance_validator import validate_provenance


def validate_evaluation_packet(
    packet: EvaluationPacket,
    sequence: str,
    protected_residues: set[int],
    entry: TargetRegistryEntry | None = None,
) -> list[str]:
    """Return deterministic validation failures; an empty list means pass."""
    failures: list[str] = []
    try:
        validate_mutation(packet.candidate.mutation, sequence, protected_residues)
    except ValueError as exc:
        failures.append(str(exc))
    if entry and packet.candidate.agi != entry.agi:
        failures.append("Candidate AGI does not match fixed target registry entry")
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
    if len(packet.scores.component_evidence) != 9:
        failures.append("Every score component must include evidence")
    if not packet.mechanistic_hypothesis:
        failures.append("Review packet is missing a mechanistic hypothesis")
    if not packet.risk_summary:
        failures.append("Review packet is missing explicit risk summaries")
    if not packet.unresolved_uncertainty:
        failures.append("Review packet must keep unresolved uncertainty visible")
    if entry and entry.agi == "AT3G62980":
        combined = " ".join(packet.known_facts + [packet.mechanistic_hypothesis]).casefold()
        if "molecular_glue" not in combined and "molecular glue" not in combined and "auxin_agonist" not in combined:
            failures.append("TIR1 review must preserve agonist/molecular-glue mechanism framing")
    return failures
