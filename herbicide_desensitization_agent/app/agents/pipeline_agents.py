from __future__ import annotations

from ..backends.interfaces import AffinityPredictionBackend, RosalindReasoningBackend
from .mutation_scoring import EvidenceAwareScoringAgent
from ..schemas.models import (
    EvaluationPacket,
    InteractionFingerprint,
    MutationCandidate,
    Pose,
    Provenance,
    ScorePacket,
    TargetProtein,
)
from ..validators.mutation_validator import validate_mutation


MOCK_PROVENANCE = [
    Provenance(
        source="mock://agent-pipeline",
        method="deterministic-placeholder",
        evidence_type="synthetic",
        notes="Not suitable for experimental decisions",
    )
]


class ConstrainedMutationAgent:
    def __init__(self, include_second_shell: bool = False) -> None:
        self.include_second_shell = include_second_shell

    ALTERNATIVES = {
        "A": ("V", "S"), "V": ("I", "A"), "I": ("L", "V"), "L": ("I", "M"), "M": ("L", "I"),
        "F": ("Y", "L"), "Y": ("F", "S"), "W": ("F", "Y"), "S": ("T", "A"), "T": ("S", "V"),
        "N": ("Q", "S"), "Q": ("N", "E"), "D": ("E", "N"), "E": ("D", "Q"), "K": ("R", "Q"),
        "R": ("K", "Q"), "H": ("N", "Q"), "C": ("S", "A"), "G": ("A", "S"), "P": ("A", "S"),
    }

    def propose(
        self,
        target: TargetProtein,
        protected: set[int],
        fingerprint: InteractionFingerprint | None = None,
    ) -> tuple[list[MutationCandidate], list[dict[str, str]]]:
        fingerprint = fingerprint or InteractionFingerprint(target.agi, [], [], [], sorted(protected), [], MOCK_PROVENANCE)
        prohibited = (protected | set(fingerprint.shared_protected)
                      | set(fingerprint.native_ligand_critical_protected)
                      | set(fingerprint.protected_by_context))
        positions = list(fingerprint.herbicide_selective_mutable)
        if self.include_second_shell:
            positions += fingerprint.second_shell_candidates
        rejected: list[dict[str, str]] = []
        candidates: list[MutationCandidate] = []
        for position in dict.fromkeys(positions):
            if not 1 <= position <= len(target.sequence):
                rejected.append({"mutation": f"position-{position}", "reason": "contact position outside target sequence"})
                continue
            source = target.sequence[position - 1]
            if position in prohibited:
                rejected.append({"mutation": f"{source}{position}?", "reason": "protected native/shared/context residue"})
                continue
            classification = (
                "HERBICIDE_SELECTIVE_MUTABLE"
                if position in fingerprint.herbicide_selective_mutable
                else "SECOND_SHELL_CANDIDATE"
            )
            for destination in self.ALTERNATIVES[source]:
                mutation = f"{source}{position}{destination}"
                validate_mutation(mutation, target.sequence, prohibited)
                candidates.append(MutationCandidate(
                    mutation=mutation,
                    agi=target.agi,
                    structure_residue=position,
                    rationale=f"Single substitution at a pose-supported {classification.lower().replace('_', ' ')} position.",
                    classification=classification,
                    provenance=fingerprint.provenance,
                    sequence_residue=position,
                    metadata={"generation_rule": "contact-guided-conservative-single-substitution"},
                ))
        if not candidates:
            rejected.append({"mutation": "none", "reason": "no permissible fingerprint-guided position available"})
        return candidates, rejected


class InteractionFingerprintAgent:
    def compare(
        self,
        target_agi: str,
        native_poses: list[Pose],
        herbicide_poses: list[Pose],
        protected: set[int],
    ) -> InteractionFingerprint:
        native_contacts = {residue for pose in native_poses for residue in pose.contacts}
        herbicide_contacts = {residue for pose in herbicide_poses for residue in pose.contacts}
        shared = native_contacts & herbicide_contacts
        herbicide_only = herbicide_contacts - native_contacts - protected
        native_only = native_contacts - herbicide_contacts
        # Sequence adjacency does not establish spatial second-shell contact.
        second_shell = {
            residue for pose in herbicide_poses
            for residue in pose.metadata.get("spatial_second_shell_residues", [])
        } - native_contacts - herbicide_contacts - protected
        return InteractionFingerprint(
            target_agi=target_agi,
            herbicide_selective_mutable=sorted(herbicide_only),
            shared_protected=sorted(shared),
            native_ligand_critical_protected=sorted(native_only),
            protected_by_context=sorted(protected),
            second_shell_candidates=sorted(second_shell),
            provenance=[
                Provenance(
                    source="computed://pose-ensemble-contact-sets",
                    method="ensemble-contact-set-comparison",
                    evidence_type="computed-fingerprint",
                )
            ],
        )


class MultiOracleScoringAgent:
    def __init__(self, affinity: AffinityPredictionBackend) -> None:
        self.affinity = affinity
        self.baseline = EvidenceAwareScoringAgent()

    def score(
        self, candidate: MutationCandidate, native: list[Pose], herbicide: list[Pose], fingerprint: InteractionFingerprint
    ) -> ScorePacket:
        packet = self.baseline.score(candidate, native, herbicide, fingerprint)
        affinity_retention = self.affinity.relative_score(native, herbicide)
        evidence = dict(packet.component_evidence)
        evidence["native_ligand_retention_score"] += f" Backend ensemble retention={affinity_retention:.3f}."
        return ScorePacket(
            **{
                field: (round((getattr(packet, field) + affinity_retention) / 2, 3) if field == "native_ligand_retention_score" else getattr(packet, field))
                for field in (
                    "herbicide_escape_score", "native_ligand_retention_score", "functional_geometry_score",
                    "fold_stability_score", "cofactor_or_complex_retention_score", "conservation_score",
                    "pose_robustness_score", "method_agreement_score", "experimental_actionability_score",
                )
            },
            provenance=packet.provenance,
            component_evidence=evidence,
            uncertainty=packet.uncertainty,
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
        uncertainty = list(dict.fromkeys(list(review["uncertainty"]) + list(scores.uncertainty)))
        return EvaluationPacket(
            candidate=candidate,
            scores=scores,
            known_facts=facts,
            predictions=predictions,
            assumptions=list(review["assumptions"]),
            unresolved_uncertainty=uncertainty,
            recommendation="more_computation",
            status="NEEDS_REVIEW",
            provenance=MOCK_PROVENANCE,
            mechanistic_hypothesis=(
                f"{candidate.mutation} may perturb a {candidate.classification.lower().replace('_', ' ')} "
                "while retaining native-function geometry; this remains a computational hypothesis."
            ),
            herbicide_interactions_disrupted=[scores.component_evidence.get("herbicide_escape_score", "")],
            native_function_interactions_preserved=[
                scores.component_evidence.get("native_ligand_retention_score", "")
            ],
            risk_summary={
                "fold": scores.component_evidence.get("fold_stability_score", "not assessed"),
                "conservation": scores.component_evidence.get("conservation_score", "not assessed"),
                "cofactor_or_complex": scores.component_evidence.get(
                    "cofactor_or_complex_retention_score", "not assessed"
                ),
            },
            recommended_assay_category="computational-validation-only",
        )
