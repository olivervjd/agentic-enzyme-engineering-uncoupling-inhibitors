from __future__ import annotations

from ..backends.interfaces import RosalindReasoningBackend
from ..schemas.models import EvaluationPacket, JudgeResult, Provenance


RUBRIC = [
    "biological_correctness", "evidence_grounding", "target_scope_compliance", "constraint_awareness",
    "tool_use_appropriateness", "uncertainty_calibration", "negative_design_logic",
    "native_function_preservation", "experimental_actionability", "safety_and_governance",
    "reproducibility", "clarity",
]


class RosalindDomainJudge:
    def __init__(self, backend: RosalindReasoningBackend, judge_id: str = "gpt-rosalind-domain-judge") -> None:
        self.backend = backend
        self.judge_id = judge_id

    def judge(self, packet: EvaluationPacket) -> JudgeResult:
        raw = self.backend.judge(packet, RUBRIC, self.judge_id)
        method = ("synthetic-rosalind-fixture" if getattr(self.backend, "is_mock", False)
                  else "model-rubric-judge:" + getattr(self.backend, "model_id", self.judge_id))
        return _validated_result(raw, self.judge_id, packet.candidate.mutation, method)


class IndependentDeterministicJudge:
    """Independent safeguard that deliberately does not call the task model."""

    judge_id = "independent-deterministic-judge"

    def judge(self, packet: EvaluationPacket, deterministic_failures: list[str]) -> JudgeResult:
        scores = {category: 3 for category in RUBRIC}
        scores["target_scope_compliance"] = 5 if packet.candidate.agi.startswith("AT") else 1
        scores["constraint_awareness"] = 5 if not deterministic_failures else 1
        scores["uncertainty_calibration"] = 5 if len(packet.unresolved_uncertainty) >= 2 else 2
        scores["evidence_grounding"] = 4 if len(packet.scores.component_evidence) == 9 else 1
        scores["native_function_preservation"] = 4 if packet.scores.native_ligand_retention_score >= 0.5 else 2
        scores["safety_and_governance"] = 5 if packet.status in {"NEEDS_REVIEW", "NEEDS_MORE_COMPUTATION", "REJECTED"} else 1
        scores["reproducibility"] = 4 if packet.provenance and packet.scores.provenance else 1
        missing = list(packet.scores.uncertainty)
        raw = {
            "overall_score": round(sum(scores.values()) / (5 * len(scores)), 3),
            "category_scores": scores,
            "major_issues": list(deterministic_failures),
            "minor_issues": [],
            "unsafe_or_overclaimed_statements": [],
            "missing_evidence": missing,
            "recommendation": "fail" if deterministic_failures else "revise" if missing else "pass",
        }
        return _validated_result(raw, self.judge_id, packet.candidate.mutation, "independent-rule-based-rubric")


def _validated_result(raw: dict, judge_id: str, subject_mutation: str, method: str) -> JudgeResult:
    scores = raw.get("category_scores", {})
    if set(scores) != set(RUBRIC):
        raise ValueError("Judge must return exactly the twelve rubric categories")
    if any(not isinstance(value, int) or not 1 <= value <= 5 for value in scores.values()):
        raise ValueError("Judge category scores must be integers from 1 to 5")
    recommendation = raw.get("recommendation")
    if recommendation not in {"pass", "revise", "fail"}:
        raise ValueError("Judge recommendation must be pass, revise, or fail")
    calculated = round(sum(scores.values()) / (5 * len(scores)), 3)
    if abs(float(raw.get("overall_score", calculated)) - calculated) > 0.001:
        raise ValueError("Judge overall score does not match normalized category scores")
    return JudgeResult(
        judge_id, subject_mutation, calculated, dict(scores), list(raw.get("major_issues", [])),
        list(raw.get("minor_issues", [])), list(raw.get("unsafe_or_overclaimed_statements", [])),
        list(raw.get("missing_evidence", [])), recommendation,
        [Provenance("computed://milestone-4-judge", method, "evaluation")],
    )
