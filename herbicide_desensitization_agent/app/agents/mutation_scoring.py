from __future__ import annotations

from dataclasses import fields

from ..schemas.models import InteractionFingerprint, MutationCandidate, ParetoRank, Pose, Provenance, ScorePacket


COMPUTED_PROVENANCE = [
    Provenance(
        source="computed://milestone-3-baseline-oracles",
        method="contact-support-and-physicochemical-heuristics",
        evidence_type="computed-baseline",
        notes="Requires replacement or corroboration by conservation, stability, and affinity tools",
    )
]

AA_PROPERTIES = {
    "A": (0, 0, 0), "V": (0, 0, 1), "I": (0, 0, 1), "L": (0, 0, 1), "M": (0, 0, 1),
    "F": (0, 1, 1), "W": (0, 1, 1), "Y": (1, 1, 1),
    "S": (1, 0, 0), "T": (1, 0, 1), "N": (1, 0, 0), "Q": (1, 0, 1), "C": (1, 0, 0),
    "D": (1, -1, 0), "E": (1, -1, 1), "K": (1, 1, 1), "R": (1, 1, 1), "H": (1, 1, 0),
    "G": (0, 0, 0), "P": (0, 0, 1),
}


def substitution_similarity(source: str, destination: str) -> float:
    left, right = AA_PROPERTIES[source], AA_PROPERTIES[destination]
    differences = sum(a != b for a, b in zip(left, right))
    return round(1.0 - differences / len(left), 3)


def contact_support(poses: list[Pose], position: int) -> float:
    if not poses:
        return 0.0
    return sum(position in pose.contacts for pose in poses) / len(poses)


class EvidenceAwareScoringAgent:
    """Transparent Milestone 3 baselines pending external physics oracles."""

    def score(
        self,
        candidate: MutationCandidate,
        native: list[Pose],
        herbicide: list[Pose],
        fingerprint: InteractionFingerprint,
    ) -> ScorePacket:
        source, destination = candidate.mutation[0], candidate.mutation[-1]
        position = candidate.sequence_residue or candidate.structure_residue
        similarity = substitution_similarity(source, destination)
        native_support = contact_support(native, position)
        herbicide_support = contact_support(herbicide, position)
        selective = position in fingerprint.herbicide_selective_mutable
        second_shell = position in fingerprint.second_shell_candidates

        escape = min(1.0, 0.20 + 0.60 * herbicide_support + (0.15 if selective else 0.05 if second_shell else 0.0))
        native_retention = max(0.0, 0.90 - native_support * (0.75 - 0.45 * similarity))
        geometry = max(0.0, 0.55 + 0.35 * similarity - 0.20 * native_support)
        predicted_ddg = candidate.metadata.get("predicted_ddg_kcal_mol")
        fold = (
            max(0.0, min(1.0, 1.0 - max(0.0, float(predicted_ddg)) / 5.0))
            if predicted_ddg is not None
            else max(0.0, 0.45 + 0.50 * similarity - (0.08 if source == "P" or destination == "P" else 0.0))
        )
        context = max(0.0, 0.85 - 0.30 * native_support) if position not in fingerprint.protected_by_context else 0.0
        homolog_frequency = candidate.metadata.get("homolog_substitution_frequency")
        conservation = float(homolog_frequency) if homolog_frequency is not None else 0.5
        robustness = max(0.0, herbicide_support - native_support / 2) if selective else 0.5 * herbicide_support
        agreement_values = [escape, native_retention, geometry, fold, context, conservation, robustness]
        agreement = 1.0 - (max(agreement_values) - min(agreement_values))
        actionability = min(native_retention, fold, context) * (0.7 + 0.3 * robustness)

        values = [escape, native_retention, geometry, fold, context, conservation, robustness, agreement, actionability]
        values = [round(max(0.0, min(1.0, value)), 3) for value in values]
        evidence = {
            "herbicide_escape_score": f"Herbicide pose contact support={herbicide_support:.2f}; selective={selective}.",
            "native_ligand_retention_score": f"Native pose contact support={native_support:.2f}; substitution similarity={similarity:.2f}.",
            "functional_geometry_score": "Proxy based on substitution similarity and native-contact support.",
            "fold_stability_score": (
                f"Converted from predicted delta-delta-G={float(predicted_ddg):.3f} kcal/mol."
                if predicted_ddg is not None else
                "Physicochemical similarity proxy; no calculated delta-delta-G is available."
            ),
            "cofactor_or_complex_retention_score": "Protected-context exclusion plus native-contact support proxy.",
            "conservation_score": (
                f"Observed homolog substitution frequency={float(homolog_frequency):.3f}."
                if homolog_frequency is not None else
                "Neutral prior because no homolog alignment was supplied."
            ),
            "pose_robustness_score": "Fraction of pose ensemble members supporting the candidate contact.",
            "method_agreement_score": "Agreement range across independent score components.",
            "experimental_actionability_score": "Minimum retention score weighted by pose robustness.",
        }
        uncertainty = ["Contact-set deltas are geometric hypotheses, not binding free energies."]
        if homolog_frequency is None:
            uncertainty.append("No homolog alignment was supplied; conservation uses a neutral prior.")
        if predicted_ddg is None:
            uncertainty.append("No Rosetta or molecular-dynamics fold calculation was supplied.")
        return ScorePacket(
            *values,
            provenance=COMPUTED_PROVENANCE,
            component_evidence=evidence,
            uncertainty=uncertainty,
        )


SCORE_FIELDS = [
    item.name for item in fields(ScorePacket)
    if item.name.endswith("_score")
]


def pareto_rank(candidates: list[MutationCandidate], scores: list[ScorePacket]) -> list[ParetoRank]:
    if len(candidates) != len(scores):
        raise ValueError("Candidates and scores must have equal lengths")
    vectors = [tuple(float(getattr(score, name)) for name in SCORE_FIELDS) for score in scores]
    remaining = set(range(len(candidates)))
    fronts: dict[int, int] = {}
    front = 1
    while remaining:
        nondominated = {
            index for index in remaining
            if not any(_dominates(vectors[other], vectors[index]) for other in remaining if other != index)
        }
        for index in nondominated:
            fronts[index] = front
        remaining -= nondominated
        front += 1
    result = []
    for index, candidate in enumerate(candidates):
        dominated_by = [candidates[j].mutation for j in range(len(candidates)) if _dominates(vectors[j], vectors[index])]
        dominates = [candidates[j].mutation for j in range(len(candidates)) if _dominates(vectors[index], vectors[j])]
        result.append(ParetoRank(candidate.mutation, fronts[index], dominated_by, dominates))
    return sorted(result, key=lambda item: (item.front, item.mutation))


def _dominates(left: tuple[float, ...], right: tuple[float, ...]) -> bool:
    return all(a >= b for a, b in zip(left, right)) and any(a > b for a, b in zip(left, right))
