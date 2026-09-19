from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from ..schemas.models import EvaluationPacket, MutationCandidate, TargetRegistryEntry, WorkflowRequest


class ScientificPlannerAgent(ABC):
    @abstractmethod
    def plan(self, request: WorkflowRequest, entry: TargetRegistryEntry) -> list[str]: ...


class EvidenceSynthesisAgent(ABC):
    @abstractmethod
    def synthesize(self, entry: TargetRegistryEntry) -> dict[str, Any]: ...


class BiologicalCriticAgent(ABC):
    @abstractmethod
    def critique(self, candidate: MutationCandidate, entry: TargetRegistryEntry) -> list[str]: ...


class CandidateReviewAgentInterface(ABC):
    @abstractmethod
    def review_candidate(self, candidate: MutationCandidate) -> EvaluationPacket: ...


class DomainJudgeAgent(ABC):
    @abstractmethod
    def judge(self, packet: EvaluationPacket) -> dict[str, Any]: ...

