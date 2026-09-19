from __future__ import annotations

from collections import defaultdict
import math

from ..registry.loader import TargetRegistry
from ..schemas.models import (
    AssayResult, EvaluationPacket, ModelRecalibrationReport, NextRoundCandidate, Provenance,
)
from ..validators.mutation_validator import MUTATION


LEARNING_PROVENANCE = [
    Provenance(
        source="computed://milestone-5-learning-loop",
        method="component-residual-recalibration-and-uncertainty-acquisition",
        evidence_type="computed-learning-output",
        notes="Requires expert review and prospective validation",
    )
]


class LearningAgent:
    def __init__(self, registry: TargetRegistry) -> None:
        self.registry = registry

    def validate_assays(
        self, assays: list[AssayResult], packets: list[EvaluationPacket], sequence: str
    ) -> None:
        candidate_mutations = {packet.candidate.mutation for packet in packets}
        seen_ids: set[str] = set()
        for assay in assays:
            entry = self.registry.get(assay.agi, assay.herbicide)
            if assay.assay_id in seen_ids:
                raise ValueError(f"Duplicate assay_id: {assay.assay_id}")
            seen_ids.add(assay.assay_id)
            match = MUTATION.fullmatch(assay.mutation)
            if not match or int(match.group(2)) > len(sequence) or sequence[int(match.group(2)) - 1] != match.group(1):
                raise ValueError(f"Invalid assay mutation for supplied sequence: {assay.mutation}")
            if assay.mutation not in candidate_mutations:
                raise ValueError(f"Assay mutation was not present in reviewed candidates: {assay.mutation}")
            required = {"herbicide_response_fold_change", "protein_expression_fraction"}
            if entry.agi == "ATCG00020":
                required.add("photosynthetic_performance_fraction")
            elif entry.agi == "AT3G62980":
                required.add("native_signaling_fraction")
            else:
                required.add("native_activity_fraction")
            missing = required - set(assay.measurements)
            if missing:
                raise ValueError(f"Assay {assay.assay_id} is missing measurements: {sorted(missing)}")
            if not assay.provenance:
                raise ValueError(f"Assay {assay.assay_id} lacks provenance")
            for name, value in assay.measurements.items():
                if not math.isfinite(float(value)):
                    raise ValueError(f"Assay measurement must be finite: {name}")
                if name.endswith("_fraction") and not 0.0 <= float(value) <= 2.0:
                    raise ValueError(f"Assay fraction outside accepted validation range [0, 2]: {name}")
                if name not in assay.units:
                    raise ValueError(f"Assay measurement lacks a unit: {name}")

    def recalibrate(
        self,
        agi: str,
        herbicide: str,
        packets: list[EvaluationPacket],
        assays: list[AssayResult],
        sequence: str,
    ) -> ModelRecalibrationReport:
        self.registry.get(agi, herbicide)
        if not assays:
            raise ValueError("At least one assay result is required for recalibration")
        if any(assay.agi != agi or assay.herbicide.casefold() != herbicide.casefold() for assay in assays):
            raise ValueError("All assays must match the recalibrated AGI-herbicide pairing")
        self.validate_assays(assays, packets, sequence)
        packet_by_mutation = {packet.candidate.mutation: packet for packet in packets}
        residuals: dict[str, list[float]] = defaultdict(list)
        absolute_errors: dict[str, list[float]] = defaultdict(list)
        for assay in assays:
            packet = packet_by_mutation[assay.mutation]
            observed = self._observed_components(assay)
            for component, value in observed.items():
                predicted = float(getattr(packet.scores, component))
                residuals[component].append(value - predicted)
                absolute_errors[component].append(abs(value - predicted))
        bias = {name: round(sum(values) / len(values), 3) for name, values in residuals.items()}
        mae = {name: round(sum(values) / len(values), 3) for name, values in absolute_errors.items()}
        calibrated: dict[str, dict[str, float]] = {}
        for packet in packets:
            calibrated[packet.candidate.mutation] = {
                component: round(_clamp(float(getattr(packet.scores, component)) + adjustment), 3)
                for component, adjustment in bias.items()
            }
        warnings = [
            "Recalibration is based on fewer than five assay records; estimates are unstable."
        ] if len(assays) < 5 else []
        warnings.append("Component residual correction is not a causal biological model.")
        return ModelRecalibrationReport(
            agi, herbicide, len(assays), bias, mae, calibrated, warnings, LEARNING_PROVENANCE
        )

    def select_next_round(
        self,
        packets: list[EvaluationPacket],
        assays: list[AssayResult],
        report: ModelRecalibrationReport,
        limit: int = 3,
    ) -> list[NextRoundCandidate]:
        if limit < 1:
            raise ValueError("Next-round candidate limit must be positive")
        assayed = {assay.mutation for assay in assays}
        ranked = []
        for packet in packets:
            mutation = packet.candidate.mutation
            if mutation in assayed:
                continue
            base = {
                "herbicide_escape_score": packet.scores.herbicide_escape_score,
                "native_ligand_retention_score": packet.scores.native_ligand_retention_score,
                "fold_stability_score": packet.scores.fold_stability_score,
                "cofactor_or_complex_retention_score": packet.scores.cofactor_or_complex_retention_score,
            }
            components = dict(base)
            components.update(report.calibrated_components.get(mutation, {}))
            exploitation = sum(components.values()) / len(components)
            exploration = min(1.0, len(packet.unresolved_uncertainty) / 4)
            acquisition = round(0.7 * exploitation + 0.3 * exploration, 3)
            reason = (
                f"Unassayed single substitution; calibrated utility={exploitation:.3f}, "
                f"uncertainty exploration={exploration:.3f}."
            )
            ranked.append(NextRoundCandidate(
                mutation, packet.candidate.agi, acquisition, reason,
                {name: round(value, 3) for name, value in components.items()},
                "NEEDS_REVIEW", LEARNING_PROVENANCE,
            ))
        return sorted(ranked, key=lambda item: (-item.acquisition_score, item.mutation))[:limit]

    @staticmethod
    def _observed_components(assay: AssayResult) -> dict[str, float]:
        measurements = assay.measurements
        native = measurements.get(
            "native_activity_fraction",
            measurements.get("native_signaling_fraction", measurements.get("photosynthetic_performance_fraction", 0.0)),
        )
        fold_change = max(0.0, measurements["herbicide_response_fold_change"])
        observed = {
            "herbicide_escape_score": fold_change / (1.0 + fold_change),
            "native_ligand_retention_score": _clamp(native),
            "cofactor_or_complex_retention_score": _clamp(measurements["protein_expression_fraction"]),
        }
        if "thermal_stability_delta_c" in measurements:
            observed["fold_stability_score"] = _clamp(0.5 + measurements["thermal_stability_delta_c"] / 10.0)
        return observed


def _clamp(value: float) -> float:
    return max(0.0, min(1.0, float(value)))
