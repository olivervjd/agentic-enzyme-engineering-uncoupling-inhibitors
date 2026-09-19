from __future__ import annotations

from ..schemas.models import StructureModel, TargetRegistryEntry


def validate_structure_context(structure: StructureModel, entry: TargetRegistryEntry) -> list[str]:
    present = {item.casefold() for item in structure.context}
    failures = [
        f"missing required context: {required}"
        for required in entry.required_context
        if required.casefold() not in present
    ]
    if entry.agi == "ATCG00020":
        for required in ("photosystem-II-complex", "thylakoid-membrane"):
            if required.casefold() not in present:
                failures.append(f"PsbA model cannot be treated as an isolated soluble protein: {required}")
    if entry.agi == "AT3G62980" and "aux-iaa-degron" not in present:
        failures.append("TIR1 model lacks Aux/IAA degron context")
    if structure.confidence < 0.0 or structure.confidence > 1.0:
        failures.append("structure confidence must be normalized to [0, 1]")
    return failures

