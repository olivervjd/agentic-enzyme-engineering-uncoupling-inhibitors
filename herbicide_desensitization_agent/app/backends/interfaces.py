from __future__ import annotations

from abc import ABC, abstractmethod

from ..schemas.models import EvaluationPacket, Ligand, MutationCandidate, Pose, StructureModel, TargetProtein, TargetRegistryEntry


class StructurePredictionBackend(ABC):
    @abstractmethod
    def predict_ensemble(self, target: TargetProtein, entry: TargetRegistryEntry) -> list[StructureModel]: ...


class DockingBackend(ABC):
    @abstractmethod
    def dock(self, structures: list[StructureModel], ligand: Ligand) -> list[Pose]: ...


class ComplexModelingBackend(ABC):
    @abstractmethod
    def validate_context(self, structure: StructureModel, entry: TargetRegistryEntry) -> list[str]: ...


class AffinityPredictionBackend(ABC):
    @abstractmethod
    def relative_score(self, native_poses: list[Pose], herbicide_poses: list[Pose]) -> float: ...


class GPUInferenceBackend(ABC):
    @abstractmethod
    def healthcheck(self) -> bool: ...


class RosalindReasoningBackend(ABC):
    @abstractmethod
    def review(self, facts: list[str], predictions: list[str]) -> dict[str, list[str] | str]: ...

    def judge(self, packet: EvaluationPacket, rubric: list[str], judge_id: str) -> dict:
        raise NotImplementedError("This Rosalind backend does not implement domain judging")


class FunctionRetentionBackend(ABC):
    """Return direct mutant structural and Kd evidence without converting IC50 to Kd."""

    def reference(
        self, target: TargetProtein, herbicide: Ligand, native_ligands: list[Ligand],
        wild_type_structures: list[StructureModel],
    ) -> dict:
        return {}

    @abstractmethod
    def evaluate(
        self,
        target: TargetProtein,
        candidate: MutationCandidate,
        herbicide: Ligand,
        native_ligands: list[Ligand],
        wild_type_structures: list[StructureModel],
    ) -> dict: ...
