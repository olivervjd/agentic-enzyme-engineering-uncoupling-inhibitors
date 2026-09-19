from __future__ import annotations

import hashlib

from ..backends.interfaces import AffinityPredictionBackend, RosalindReasoningBackend
from ..schemas.models import EvaluationPacket, MutationCandidate, Pose, Provenance, ScorePacket, TargetProtein
from ..validators.mutation_validator import validate_mutation


MOCK_PROVENANCE = [
    Provenance(
        source="mock://agent-pipeline",
        method="deterministic-placeholder",
        evidence_type="synthetic",
        notes="Not suitable for experimental decisions",
    )
]


def _score(label: str, low: float = 0.35, high: float = 0.75) -> float:
    fraction = int(hashlib.sha256(label.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
    return round(low + fraction * (high - low), 3)


class ConstrainedMutationAgent:
    def propose(self, target: TargetProtein, protected: set[int]) -> tuple[list[MutationCandidate], list[dict[str, str]]]:
        alternatives = {"A": "V", "V": "I", "I": "L", "L": "I", "M": "L"}
        rejected: list[dict[str, str]] = []
        for position, source in enumerate(target.sequence, 1):
            if position in protected:
                rejected.append({"mutation": f"{source}{position}?", "reason": "protected residue"})
                continue
            destination = alternatives.get(source, "A" if source != "A" else "V")
            mutation = f"{source}{position}{destination}"
            validate_mutation(mutation, target.sequence, protected)
            return [
                MutationCandidate(
                    mutation=mutation,
                    agi=target.agi,
                    structure_residue=position,
                    rationale="Synthetic conservative substitution used to exercise workflow plumbing only.",
                    classification="SECOND_SHELL_CANDIDATE",
                    provenance=MOCK_PROVENANCE,
                )
            ], rejected
        return [], rejected + [{"mutation": "none", "reason": "no unprotected residue available"}]


class MultiOracleScoringAgent:
    def __init__(self, affinity: AffinityPredictionBackend) -> None:
        self.affinity = affinity

    def score(self, candidate: MutationCandidate, native: list[Pose], herbicide: list[Pose]) -> ScorePacket:
        label = candidate.agi + candidate.mutation
        return ScorePacket(
            herbicide_escape_score=_score(label + "escape"),
            native_ligand_retention_score=self.affinity.relative_score(native, herbicide),
            functional_geometry_score=_score(label + "geometry"),
            fold_stability_score=_score(label + "fold"),
            cofactor_or_complex_retention_score=_score(label + "context"),
            conservation_score=_score(label + "conservation"),
            pose_robustness_score=_score(label + "pose"),
            method_agreement_score=_score(label + "agreement"),
            experimental_actionability_score=_score(label + "actionability"),
            provenance=MOCK_PROVENANCE,
        )


class CandidateReviewAgent:
    def __init__(self, reasoner: RosalindReasoningBackend) -> None:
        self.reasoner = reasoner

    def review(self, candidate: MutationCandidate, scores: ScorePacket, facts: list[str]) -> EvaluationPacket:
        predictions = [
            f"Synthetic score packet generated for {candidate.mutation}.",
            "No score represents a validated biological prediction.",
        ]
        review = self.reasoner.review(facts, predictions)
        return EvaluationPacket(
            candidate=candidate,
            scores=scores,
            known_facts=facts,
            predictions=predictions,
            assumptions=list(review["assumptions"]),
            unresolved_uncertainty=list(review["uncertainty"]),
            recommendation="more_computation",
            status="NEEDS_REVIEW",
            provenance=MOCK_PROVENANCE,
        )

