from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from .interfaces import (
    AffinityPredictionBackend, ComplexModelingBackend, DockingBackend, FunctionRetentionBackend,
    StructurePredictionBackend,
)
from ..schemas.models import Ligand, MutationCandidate, Pose, Provenance, StructureModel, TargetProtein, TargetRegistryEntry


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


class PrecomputedFunctionRetentionBackend(FunctionRetentionBackend):
    """Replay direct Kd and mutant-structure measurements with explicit provenance."""

    def __init__(self, records: dict[str, dict], base_directory: str | Path = ".") -> None:
        self.records = records
        self.base_directory = Path(base_directory).resolve()

    def reference(self, target, herbicide, native_ligands, wild_type_structures):
        return self._record("WT", target)

    def evaluate(
        self, target: TargetProtein, candidate: MutationCandidate, herbicide: Ligand,
        native_ligands: list[Ligand], wild_type_structures: list[StructureModel],
    ) -> dict:
        return self._record(candidate.mutation, target)

    def _record(self, mutation, target):
        record = deepcopy(self.records.get(mutation, {}))
        if record and record.get("target_agi") != target.agi:
            raise ValueError("Retention evidence requires the matching target_agi")
        provenance = record.get("provenance", [])
        record["provenance"] = [
            item if isinstance(item, Provenance) else Provenance(**item) for item in provenance
        ]
        comparisons = record.get("structure_comparisons", [])
        if comparisons:
            from .structural import TMAlignStructuralMatcher

            matcher = TMAlignStructuralMatcher()
            metrics = []
            for comparison in comparisons:
                options = dict(comparison)
                for name in ("mutant", "reference"):
                    options[name] = self.base_directory / options[name]
                metrics.append(matcher.compare(**options, target_sequence=target.sequence, mutation=mutation))
                record["provenance"].append(Provenance(
                    str(options["mutant"]), matcher.method, "computed-structure-comparison",
                    f"Reference: {options['reference']}; parameters: {comparison}",
                ))
            # Every supplied condition/replicate must pass; do not cherry-pick the best model.
            record["structural_metrics"] = {
                key: None if any(item[key] is None for item in metrics) else
                (max(item[key] for item in metrics) if "rmsd" in key else min(item[key] for item in metrics))
                for key in metrics[0]
            }
            record["structural_method"] = matcher.method + "; worst case over supplied comparisons"
        return record
