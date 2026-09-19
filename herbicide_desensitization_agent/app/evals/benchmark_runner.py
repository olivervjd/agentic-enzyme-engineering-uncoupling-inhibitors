from __future__ import annotations

from dataclasses import dataclass

from ..registry.loader import TargetRegistry
from ..schemas.models import EvaluationPacket
from ..validators.mutation_validator import validate_mutation


@dataclass(frozen=True)
class BenchmarkCase:
    case_id: str
    agi: str
    herbicide: str
    sequence: str
    protected_residues: set[int]
    withheld_known_positions: set[int]
    evidence_type: str


class FixedTargetBenchmarkRunner:
    def __init__(self, registry: TargetRegistry) -> None:
        self.registry = registry

    def run(self, case: BenchmarkCase, packets: list[EvaluationPacket]) -> dict:
        self.registry.get(case.agi, case.herbicide)
        if case.evidence_type not in {"synthetic", "curated-retrospective"}:
            raise ValueError("Benchmark labels must be disclosed as synthetic or curated-retrospective")
        ranked = [packet.candidate.sequence_residue or packet.candidate.structure_residue for packet in packets]
        known = case.withheld_known_positions
        invalid = 0
        protected = 0
        for packet in packets:
            try:
                _, position, _ = validate_mutation(packet.candidate.mutation, case.sequence, set())
                protected += int(position in case.protected_residues)
            except ValueError:
                invalid += 1
        count = max(1, len(packets))
        denominator = max(1, len(known))
        return {
            "case_id": case.case_id, "agi": case.agi, "herbicide": case.herbicide,
            "evidence_type": case.evidence_type,
            "top_5_known_position_recall": round(len(set(ranked[:5]) & known) / denominator, 3),
            "top_10_known_position_recall": round(len(set(ranked[:10]) & known) / denominator, 3),
            "protected_residue_violation_rate": round(protected / count, 3),
            "invalid_mutation_rate": round(invalid / count, 3),
            "native_function_retention_quality": _mean(p.scores.native_ligand_retention_score for p in packets),
            "evidence_completeness_score": _mean(len(p.scores.component_evidence) / 9 for p in packets),
            "uncertainty_quality_score": _mean(min(1.0, len(p.unresolved_uncertainty) / 2) for p in packets),
            "assay_actionability_score": _mean(p.scores.experimental_actionability_score for p in packets),
        }


def _mean(values) -> float:
    values = list(values)
    return round(sum(values) / len(values), 3) if values else 0.0
