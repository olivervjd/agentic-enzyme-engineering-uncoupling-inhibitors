from __future__ import annotations

from dataclasses import asdict
from typing import Callable

from .interfaces import ReasoningBackend
from ..agents.prompts import CANDIDATE_REVIEW_PROMPT, DOMAIN_JUDGE_PROMPT
from ..schemas.models import EvaluationPacket


class GPTJSONBackend(ReasoningBackend):
    """Transport-injected structured reasoning adapter with no persisted credentials."""

    def __init__(self, transport: Callable[[str, dict], dict]) -> None:
        self.transport = transport

    def review(self, facts: list[str], predictions: list[str]) -> dict:
        result = self.transport("candidate_review", {
            "system_prompt": CANDIDATE_REVIEW_PROMPT,
            "facts": facts,
            "predictions": predictions,
            "required_keys": ["assumptions", "uncertainty", "recommendation"],
        })
        for key in ("assumptions", "uncertainty", "recommendation"):
            if key not in result:
                raise ValueError(f"Model review response is missing {key}")
        return result

    def judge(self, packet: EvaluationPacket, rubric: list[str], judge_id: str) -> dict:
        return self.transport("domain_judge", {
            "system_prompt": DOMAIN_JUDGE_PROMPT,
            "judge_id": judge_id,
            "rubric": rubric,
            "packet": asdict(packet),
            "required_recommendations": ["pass", "revise", "fail"],
        })
