"""Fail-closed calibration and binding-site gates for live workflows."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


class WorkflowReassessmentRequired(ValueError):
    def __init__(self, reasons):
        self.reasons = list(reasons)
        super().__init__("Workflow needs reassessment: " + "; ".join(self.reasons))


class StructureQualityAgent:
    def __init__(self, report_path, binding_report_path=None):
        self.path = Path(report_path)
        self.binding_path = Path(binding_report_path) if binding_report_path else None

    def _load(self, request):
        report = json.loads(self.path.read_text())
        expected = hashlib.sha256(request.target.sequence.encode()).hexdigest()
        if report.get("target_agi") != request.target.agi or report.get("sequence_sha256") != expected:
            raise ValueError("Calibration target/sequence does not match workflow request")
        for artifact in report.get("artifacts", []):
            path = self.path.parent / artifact["path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
                raise ValueError("Calibration artifact changed: " + artifact["path"])
        return report

    def assess(self, request, structures):
        report = self._load(request)
        checks = report.get("checks", {})
        required = ("msa_validated", "chemistry_reviewed", "wt_repeatability", "experimental_accuracy",
                    "functional_controls", "no_forced_templates")
        reasons = [name + " not established" for name in required if checks.get(name) is not True]
        paths = {str((self.path.parent / a["path"]).resolve()) for a in report.get("artifacts", [])}
        if not paths or any(not s.artifact_path or str(Path(s.artifact_path).resolve()) not in paths for s in structures):
            reasons.append("Structure ensemble is not bound to the calibration artifacts")
        reasons.extend(report.get("blocking_reasons", []))
        binding = None
        if self.binding_path:
            binding = json.loads(self.binding_path.read_text())
            if (binding.get("target_id") != request.target.agi or binding.get("prediction_sha256") !=
                    hashlib.sha256((self.path.parent / "predictions.json").read_bytes()).hexdigest()):
                raise ValueError("Binding calibration target/prediction mismatch")
            if binding.get("status") != "VALIDATED_WITHIN_SCOPE":
                reasons.append("Endpoint-matched binding calibration not validated")
        return {"status": "NEEDS_REASSESSMENT" if reasons else "PASSED", "reasons": sorted(set(reasons)),
                "report": str(self.path), "binding_calibration": binding}

    def assess_pocket(self, request, fingerprint, native, herbicide, diagnostics=None):
        report = self._load(request)
        pocket = report.get("pocket", {})
        reasons = []
        if pocket.get("matched_context") is not True:
            reasons.append("Independent pose methods do not have a validated matching chemical context")
        if diagnostics is not None and diagnostics.get("matched_context") is not True:
            reasons.append("New docking diagnostics omit required cosubstrate context")
        supported = set(pocket.get("supported_herbicide_positions", []))
        disputed = set(fingerprint.herbicide_selective_mutable) - supported
        if disputed:
            reasons.append("Uncorroborated mutation sites: " + ", ".join(map(str, sorted(disputed))))
        if not supported:
            reasons.append("No independently supported herbicide pocket")
        return {"status": "NEEDS_REASSESSMENT" if reasons else "PASSED", "reasons": reasons}
