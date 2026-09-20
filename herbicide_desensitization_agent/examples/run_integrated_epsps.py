"""Run the original orchestrator with live calibration, evidence and independent model settings."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from ..app.agents.evidence import CitedEvidenceAgent
from ..app.agents.quality import StructureQualityAgent, WorkflowReassessmentRequired
from ..app.agents.pose_ensemble import PoseEnsembleAgent
from ..app.agents.workflow_review import WorkflowReviewAgent
from ..app.backends.calibration import EPSPSCalibrationRunner, write
from ..app.backends.calibrated_epsps import ArtifactCandidateOracles, CalibratedEPSPSBackend
from ..app.backends.openai_json import DEFAULT_EVIDENCE_MODEL, configured_reasoner, OpenAIJSONTransport, UnavailableModelTransport
from ..app.backends.literature import EuropePMCRetriever
from ..app.backends.credentials import resolve_api_key
from ..app.orchestrator.workflow import WorkflowOrchestrator
from ..app.registry.loader import TargetRegistry
from ..app.schemas.models import AssayResult, FunctionalPartner, Ligand, Provenance, TargetProtein, WorkflowRequest
from ..app.storage.artifact_store import ArtifactStore
from .check_model_access import inspect_model_access


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boltz-cache", type=Path, required=True)
    parser.add_argument("--boltz", default="boltz")
    parser.add_argument("--diffdock-cache", type=Path)
    evidence_options = parser.add_mutually_exclusive_group()
    evidence_options.add_argument("--evidence-model", default=os.getenv("EVIDENCE_MODEL", DEFAULT_EVIDENCE_MODEL))
    evidence_options.add_argument("--curated-evidence-only", action="store_true",
                                  help="Explicitly disable model synthesis; never an automatic fallback")
    parser.add_argument("--review-model", default=os.getenv("REVIEW_MODEL", "gpt-5.6-luna"))
    parser.add_argument("--judge-model", default=os.getenv("JUDGE_MODEL", "gpt-5.6-luna"))
    parser.add_argument("--use-codex-api-key", action="store_true", help="Reuse only an actual locally stored Codex API key, never OAuth")
    parser.add_argument("--binding-dataset", type=Path, help="Curated endpoint-matched experimental calibration dataset")
    parser.add_argument("--diagnostics-only", action="store_true", help="Review saved predictions; never generate mutations or launch Boltz")
    parser.add_argument("--cached-only", action="store_true", help="Fail if predictions are missing; never launch Boltz")
    parser.add_argument("--no-literature-retrieval", action="store_true", help="Explicitly disable Europe PMC abstract retrieval")
    parser.add_argument("--assays", type=Path, help="JSON list of real AssayResult records with provenance")
    args = parser.parse_args()
    if args.curated_evidence_only:
        args.evidence_model = None
    if args.output.exists():
        parser.error("Use a fresh workflow output directory")
    cached_only = args.cached_only or args.diagnostics_only
    if cached_only and not all((args.calibration / name).is_file() for name in ("calibration_inputs.json", "predictions.json", "calibration_report.json")):
        parser.error("Cached-only execution requires existing calibration inputs, predictions and report")
    args.output.mkdir(parents=True)
    runner = EPSPSCalibrationRunner(args.previous, args.calibration, args.boltz_cache, args.boltz)
    api_key, credential_source = resolve_api_key(use_codex=args.use_codex_api_key)
    model_settings = {role: {"model": model, "configured": bool(model and api_key), "credential_source": credential_source}
                      for role, model in (("evidence", args.evidence_model), ("review", args.review_model), ("judge", args.judge_model))}
    access_cache = {}
    for setting in model_settings.values():
        model = setting["model"]
        if model and api_key:
            if model not in access_cache:
                try:
                    access_cache[model] = inspect_model_access(model, api_key)
                except OSError:
                    access_cache[model] = {"metadata_accessible": False, "status": "CHECK_FAILED"}
            setting.update(access_cache[model])
        else:
            setting.update(metadata_accessible=False, status="DISABLED" if not model else "NO_CREDENTIAL")
    write(args.output / "model_configuration.json", model_settings)

    def transport(role):
        setting = model_settings[role]
        if not setting["metadata_accessible"]:
            return UnavailableModelTransport(setting["model"], role + " model access unavailable; see model_configuration.json")
        return OpenAIJSONTransport(setting["model"], api_key=api_key)

    try:
        evidence_transport = None if args.curated_evidence_only else transport("evidence")
        manifest = runner.prepare() if not args.calibration.exists() else json.loads((args.calibration / "calibration_inputs.json").read_text())
        registry = TargetRegistry.default()
        entry = registry.get("AT2G45300", "glyphosate")
        prov = [Provenance(str(args.calibration / "calibration_inputs.json"), "validated-input-manifest", "reference-input")]
        request = WorkflowRequest(
            TargetProtein(entry.agi, entry.protein_name, manifest["sequence"], provenance=prov),
            Ligand("glyphosate", "herbicide", manifest["ligands"]["glyphosate"], provenance=prov),
            [Ligand("phosphoenolpyruvate", "native", manifest["ligands"]["pep"], provenance=prov),
             Ligand("shikimate-3-phosphate", "native", manifest["ligands"]["s3p"], provenance=prov)],
            [FunctionalPartner(name, "required-but-not-automatically-satisfied", prov) for name in entry.required_context],
        )
        binding_path = args.output / "binding_calibration.json"
        backend = CalibratedEPSPSBackend(runner, binding_dataset=args.binding_dataset, binding_output=binding_path,
                                         cached_only=cached_only)
        review = configured_reasoner(args.review_model, api_key=api_key) if model_settings["review"]["metadata_accessible"] else None
        judge = configured_reasoner(args.judge_model, api_key=api_key) if model_settings["judge"]["metadata_accessible"] else None
        assays = []
        if args.assays:
            for row in json.loads(args.assays.read_text()):
                assays.append(AssayResult(**{**row, "provenance": [Provenance(**p) for p in row["provenance"]]}))
        oracles = ArtifactCandidateOracles(args.calibration / "candidate_oracles.json",
                                          args.calibration / "calibration_inputs.json")
        workflow = WorkflowOrchestrator(
            registry, backend, backend, backend, backend, review, ArtifactStore(args.output),
            function_retention_backend=oracles,
            judge_backend=judge, evidence_agent=CitedEvidenceAgent(args.calibration / "evidence.json", evidence_transport,
                None if args.no_literature_retrieval else EuropePMCRetriever(args.output / "literature")),
            quality_agent=StructureQualityAgent(args.calibration / "calibration_report.json", binding_path),
            candidate_evidence_backend=oracles,
            pose_ensemble_agent=PoseEnsembleAgent(args.output / "diffdock", args.diffdock_cache) if args.diffdock_cache else None,
            require_live_gates=True,
            diagnostic_review_agent=WorkflowReviewAgent(transport("review"), transport("judge")),
            diagnostics_only=args.diagnostics_only, assays=assays,
        )
        workflow.run(request)
        incomplete = [s["stage"] for s in workflow.stage_manifest if s["status"] in {"FAILED", "UNAVAILABLE", "BLOCKED", "NOT_RUN"}]
        write(args.output / "execution.json", {"status": "COMPUTATIONAL_STAGES_COMPLETE" if not incomplete
                                               else "PARTIAL_AGENT_EXECUTION", "incomplete_stages": incomplete,
                                               "limitations": ["No automatic assay approval", "Learning awaits real assays"],
                                               "models": model_settings})
    except WorkflowReassessmentRequired as exc:
        write(args.output / "execution.json", {"status": "NEEDS_REASSESSMENT", "reasons": exc.reasons,
                                               "models": model_settings, "candidates_advanced": 0})
        print(str(exc))
        return 2
    except Exception as exc:
        write(args.output / "execution.json", {"status": "FAILED", "error_type": type(exc).__name__, "models": model_settings})
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
