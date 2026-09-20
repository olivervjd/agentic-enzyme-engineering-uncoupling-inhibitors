"""Live calibration outputs adapted to the original WorkflowOrchestrator contracts."""
from __future__ import annotations

import json
import hashlib
from dataclasses import replace
from pathlib import Path

from .interfaces import AffinityPredictionBackend, ComplexModelingBackend, DockingBackend, StructurePredictionBackend, FunctionRetentionBackend
from .precomputed import PrecomputedFunctionRetentionBackend
from .calibration import sha, write
from .binding_calibration import build_binding_report
from .structural import TMAlignStructuralMatcher
from ..agents.quality import WorkflowReassessmentRequired
from ..schemas.models import Pose, Provenance, StructureModel


LIGAND_NAMES = {"pep": "phosphoenolpyruvate", "s3p": "shikimate-3-phosphate", "glyphosate": "glyphosate"}


class CalibratedEPSPSBackend(StructurePredictionBackend, DockingBackend, ComplexModelingBackend, AffinityPredictionBackend):
    def __init__(self, runner, *, binding_dataset=None, binding_output=None, cached_only=False):
        self.runner = runner
        self.root = runner.root
        self._predictions = []
        self.binding_dataset, self.binding_output = binding_dataset, binding_output
        self.cached_only = cached_only

    def predict_ensemble(self, target, entry):
        if target.agi != "AT2G45300" or entry.agi != target.agi:
            raise ValueError("This live backend supports EPSPS only")
        if not (self.root / "predictions.json").exists():
            if self.cached_only:
                raise ValueError("Cached-only execution requires existing predictions; no GPU job started")
            self.runner.run()
        manifest = json.loads((self.root / "calibration_inputs.json").read_text())
        if manifest["sequence"] != target.sequence:
            raise ValueError("Calibration sequence mismatch")
        self._predictions = json.loads((self.root / "predictions.json").read_text())["records"]
        self._verify_artifacts()
        if self.binding_output:
            write(self.binding_output, build_binding_report(self.root, self.binding_dataset))
        models = []
        for row in self._predictions:
            if row["subject"] != "arabidopsis_wt":
                continue
            path = self.root / row["structure"]
            confidence_path = path.with_name(f"confidence_{row['id']}_model_0.json")
            confidence = json.loads(confidence_path.read_text())
            models.append(StructureModel(
                row["id"], target.agi, "Boltz-2 2.2.1 with validated MSA", confidence["confidence_score"],
                ["shikimate-3-phosphate"],
                [Provenance(str(path), "upstream-boltz2-msa-calibration", "predicted-structure",
                            "S3P modeled; ligand protonation, buffer conditions and catalytic waters require explicit validation. EPSPS has no mandatory catalytic-metal requirement in the curated baseline.")],
                str(path), "mmcif", {"seed": row["seed"], "context": row["context"], "confidence": confidence},
            ))
        return models

    def _verify_artifacts(self):
        report = json.loads((self.root / "calibration_report.json").read_text())
        for item in report["artifacts"]:
            if sha(self.root / item["path"]) != item["sha256"]:
                raise ValueError("Calibration artifact integrity failure")

    def dock(self, structures, ligand):
        ids = {s.model_id for s in structures}
        poses = []
        for row in self._predictions:
            if row["id"] not in ids:
                continue
            primary = "glyphosate" if row["context"] == "herbicide" else "pep"
            for chain, name in (("B", primary), ("C", "s3p")):
                if LIGAND_NAMES[name] != ligand.name:
                    continue
                poses.append(Pose(
                    row["id"] + "-" + chain, row["id"], ligand.name,
                    next(s.confidence for s in structures if s.model_id == row["id"]),
                    [p + 76 for p in row["contacts"][chain]],
                    [Provenance(str(self.root / row["structure"]), "cofold-heavy-atom-contact-5A", "predicted-pose")],
                    str(self.root / row["structure"]), 1,
                    {"method": "Boltz cofold, not independent docking", "context": row["context"], "seed": row["seed"]},
                ))
        return poses

    def validate_context(self, structure, entry):
        return ["missing:" + name for name in entry.required_context if name not in structure.context]

    def relative_score(self, native_poses, herbicide_poses):
        return None


class ArtifactCandidateOracles(FunctionRetentionBackend):
    """Import only protocol-matched mutant measurements, never WT proxies for mutants."""

    def __init__(self, manifest_path, protocol_sha256):
        self.path, self.protocol = Path(manifest_path), protocol_sha256
        self.retention = None

    def evaluate_candidates(self, request, candidates):
        if not self.path.exists():
            raise WorkflowReassessmentRequired(["No protocol-matched mutant oracle artifacts supplied"])
        data = json.loads(self.path.read_text())
        protocol = sha(self.protocol) if isinstance(self.protocol, Path) else self.protocol
        if data.get("target_agi") != request.target.agi or data.get("protocol_sha256") != protocol:
            raise ValueError("Candidate oracle target/protocol mismatch")
        rows = data.get("candidates", {})
        results = []
        allowed_artifacts = set()
        for candidate in candidates:
            if candidate.mutation not in rows:
                raise WorkflowReassessmentRequired(["Missing mutant oracles: " + candidate.mutation])
            row = rows[candidate.mutation]
            if not row.get("artifacts") or not row.get("measurements"):
                raise ValueError("Oracle measurements require artifact provenance")
            sequence = request.target.sequence
            position = int(candidate.mutation[1:-1])
            expected = sequence[:position - 1] + candidate.mutation[-1] + sequence[position:]
            if row.get("sequence_sha256") != hashlib.sha256(expected.encode()).hexdigest():
                raise ValueError("Mutant oracle sequence identity mismatch")
            structures = []
            for artifact in row["artifacts"]:
                if sha(self.path.parent / artifact["path"]) != artifact["sha256"]:
                    raise ValueError("Candidate oracle artifact changed")
                allowed_artifacts.add(artifact["path"])
                if artifact.get("role") == "mutant_structure":
                    residues = TMAlignStructuralMatcher._residues(
                        self.path.parent / artifact["path"], artifact["chain"], artifact["residue_offset"])
                    TMAlignStructuralMatcher._validate_sequence(
                        residues, sequence, candidate.mutation, artifact["sequence_start"], artifact["sequence_end"])
                    structures.append(artifact)
            if not structures:
                raise ValueError("Mutant oracle requires a sequence-validated mutant structure")
            results.append(replace(candidate, metadata={**candidate.metadata, "oracle_evidence": row["measurements"]}))
        retention = data.get("retention_records", {})
        for record in retention.values():
            for comparison in record.get("structure_comparisons", []):
                if not {comparison["mutant"], comparison["reference"]} <= allowed_artifacts:
                    raise ValueError("Retention comparison uses artifacts outside the validated oracle manifest")
        self.retention = PrecomputedFunctionRetentionBackend(retention, self.path.parent)
        return results

    def reference(self, target, herbicide, native_ligands, wild_type_structures):
        return self.retention.reference(target, herbicide, native_ligands, wild_type_structures) if self.retention else {}

    def evaluate(self, target, candidate, herbicide, native_ligands, wild_type_structures):
        return self.retention.evaluate(target, candidate, herbicide, native_ligands, wild_type_structures) if self.retention else {}
