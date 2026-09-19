from __future__ import annotations

import hashlib

from .interfaces import (
    AffinityPredictionBackend,
    ComplexModelingBackend,
    DockingBackend,
    GPUInferenceBackend,
    RosalindReasoningBackend,
    StructurePredictionBackend,
)
from ..schemas.models import Ligand, Pose, Provenance, StructureModel, TargetProtein, TargetRegistryEntry


def _fraction(label: str) -> float:
    return int(hashlib.sha256(label.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF


class MockScientificBackend(
    StructurePredictionBackend,
    DockingBackend,
    ComplexModelingBackend,
    AffinityPredictionBackend,
    GPUInferenceBackend,
    RosalindReasoningBackend,
):
    """Deterministic synthetic backend used only to exercise orchestration."""

    is_mock = True

    provenance = [
        Provenance(
            source="mock://scientific-backend",
            method="deterministic-placeholder",
            evidence_type="synthetic",
            notes="Not biological evidence",
        )
    ]

    def predict_ensemble(self, target: TargetProtein, entry: TargetRegistryEntry) -> list[StructureModel]:
        return [
            StructureModel(
                model_id=f"{target.agi}-mock-{index}",
                target_agi=target.agi,
                method="mock-structure",
                confidence=round(0.45 + 0.05 * index, 3),
                context=list(entry.required_context),
                provenance=self.provenance,
            )
            for index in (1, 2)
        ]

    def dock(self, structures: list[StructureModel], ligand: Ligand) -> list[Pose]:
        poses = []
        for index, structure in enumerate(structures, 1):
            seed = _fraction(structure.model_id + ligand.name)
            contacts = sorted({1 + int(seed * 4), 5 + index})
            poses.append(
                Pose(
                    pose_id=f"{structure.model_id}-{ligand.role}-{index}",
                    model_id=structure.model_id,
                    ligand_name=ligand.name,
                    confidence=round(0.35 + 0.4 * seed, 3),
                    contacts=contacts,
                    provenance=self.provenance,
                )
            )
        return poses

    def validate_context(self, structure: StructureModel, entry: TargetRegistryEntry) -> list[str]:
        missing = sorted(set(entry.required_context) - set(structure.context))
        return [f"missing:{item}" for item in missing]

    def relative_score(self, native_poses: list[Pose], herbicide_poses: list[Pose]) -> float:
        native = sum(p.confidence for p in native_poses) / len(native_poses)
        herbicide = sum(p.confidence for p in herbicide_poses) / len(herbicide_poses)
        return round(max(0.0, min(1.0, 0.5 + native - herbicide)), 3)

    def healthcheck(self) -> bool:
        return True

    def review(self, facts: list[str], predictions: list[str]) -> dict[str, list[str] | str]:
        return {
            "assumptions": ["All structure, pose, mutation, and score outputs are synthetic placeholders."],
            "uncertainty": ["No experimental or real computational evidence has been supplied."],
            "recommendation": "more_computation",
        }

    def judge(self, packet, rubric: list[str], judge_id: str) -> dict:
        scores = {category: 3 for category in rubric}
        scores["uncertainty_calibration"] = 5 if packet.unresolved_uncertainty else 1
        scores["safety_and_governance"] = 5 if packet.status == "NEEDS_REVIEW" else 1
        return {
            "overall_score": round(sum(scores.values()) / (5 * len(scores)), 3),
            "category_scores": scores,
            "major_issues": [],
            "minor_issues": ["Judgment uses a deterministic synthetic Rosalind fixture."],
            "unsafe_or_overclaimed_statements": [],
            "missing_evidence": list(packet.unresolved_uncertainty),
            "recommendation": "revise",
        }
