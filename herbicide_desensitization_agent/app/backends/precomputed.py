from __future__ import annotations

from .interfaces import AffinityPredictionBackend, ComplexModelingBackend, DockingBackend, StructurePredictionBackend
from ..schemas.models import Ligand, Pose, StructureModel, TargetProtein, TargetRegistryEntry


class PrecomputedComplexBackend(
    StructurePredictionBackend, DockingBackend, ComplexModelingBackend, AffinityPredictionBackend
):
    """Replay provenance-bearing structure and pose ensembles from a completed live run."""

    def __init__(
        self,
        target_agi: str,
        structures: list[StructureModel],
        poses_by_ligand: dict[str, list[Pose]],
    ) -> None:
        if len(structures) < 2:
            raise ValueError("A precomputed structure ensemble requires at least two models")
        self.target_agi = target_agi
        self.structures = structures
        self.poses_by_ligand = {name.casefold(): poses for name, poses in poses_by_ligand.items()}

    def predict_ensemble(self, target: TargetProtein, entry: TargetRegistryEntry) -> list[StructureModel]:
        if target.agi != self.target_agi or entry.agi != self.target_agi:
            raise ValueError("Precomputed backend target does not match request")
        return list(self.structures)

    def dock(self, structures: list[StructureModel], ligand: Ligand) -> list[Pose]:
        try:
            poses = self.poses_by_ligand[ligand.name.casefold()]
        except KeyError as exc:
            raise ValueError(f"No precomputed pose ensemble for {ligand.name}") from exc
        valid_models = {structure.model_id for structure in structures}
        return [pose for pose in poses if pose.model_id in valid_models]

    def validate_context(self, structure: StructureModel, entry: TargetRegistryEntry) -> list[str]:
        missing = sorted(set(entry.required_context) - set(structure.context))
        return [f"missing:{item}" for item in missing]

    def relative_score(self, native_poses: list[Pose], herbicide_poses: list[Pose]) -> float:
        if not native_poses or not herbicide_poses:
            return 0.0
        native = sum(pose.confidence for pose in native_poses) / len(native_poses)
        herbicide = sum(pose.confidence for pose in herbicide_poses) / len(herbicide_poses)
        return round(max(0.0, min(1.0, 0.5 + native - herbicide)), 3)
