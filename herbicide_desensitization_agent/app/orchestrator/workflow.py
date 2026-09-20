from __future__ import annotations

import json
from dataclasses import asdict, replace

from ..agents.pipeline_agents import (
    CandidateReviewAgent,
    ConstrainedMutationAgent,
    InteractionFingerprintAgent,
    MultiOracleScoringAgent,
)
from ..agents.mutation_scoring import pareto_rank
from ..agents.function_retention import FunctionRetentionAgent, RETENTION_LEGEND, retention_csv, retention_markdown
from ..backends.interfaces import (
    AffinityPredictionBackend,
    ComplexModelingBackend,
    DockingBackend,
    FunctionRetentionBackend,
    ReasoningBackend,
    StructurePredictionBackend,
)
from ..evals.deterministic_checks import validate_evaluation_packet
from ..evals.judges import IndependentDeterministicJudge, DomainJudge
from ..registry.loader import TargetRegistry
from ..schemas.models import WorkflowRequest, WorkflowResult
from ..storage.artifact_store import ArtifactStore
from ..validators.input_validator import validate_request
from ..validators.provenance_validator import validate_provenance
from ..validators.structure_context_validator import validate_structure_context
from ..visualization import MolstarArtifactRenderer


class WorkflowOrchestrator:
    def __init__(
        self,
        registry: TargetRegistry,
        structure_backend: StructurePredictionBackend,
        docking_backend: DockingBackend,
        complex_backend: ComplexModelingBackend,
        affinity_backend: AffinityPredictionBackend,
        reasoning_backend: ReasoningBackend | None,
        artifact_store: ArtifactStore | None = None,
        function_retention_backend: FunctionRetentionBackend | None = None,
        *,
        judge_backend: ReasoningBackend | None = None,
        evidence_agent=None,
        quality_agent=None,
        candidate_evidence_backend=None,
        pose_ensemble_agent=None,
        require_live_gates: bool = False,
        diagnostic_review_agent=None,
        diagnostics_only: bool = False,
        assays=None,
        evidence_calibration=None,
    ) -> None:
        self.registry = registry
        self.structure_backend = structure_backend
        self.docking_backend = docking_backend
        self.complex_backend = complex_backend
        self.mutation_agent = ConstrainedMutationAgent()
        self.fingerprint_agent = InteractionFingerprintAgent()
        self.scoring_agent = MultiOracleScoringAgent(affinity_backend)
        self.review_agent = CandidateReviewAgent(reasoning_backend)
        self.evidence_agent = evidence_agent
        self.quality_agent = quality_agent
        self.candidate_evidence_backend = candidate_evidence_backend
        self.pose_ensemble_agent = pose_ensemble_agent
        self.require_live_gates = require_live_gates
        self.diagnostic_review_agent = diagnostic_review_agent
        self.diagnostics_only = diagnostics_only
        self.assays = assays or []
        self.evidence_calibration = evidence_calibration
        self.stage_manifest = []
        selected_judge = judge_backend if require_live_gates else (judge_backend if judge_backend is not None else reasoning_backend)
        judge_id = (
            "mock-domain-judge"
            if getattr(selected_judge, "is_mock", False)
            else getattr(selected_judge, "model_id", "configured-domain-judge")
        )
        self.domain_judge = DomainJudge(selected_judge, judge_id) if selected_judge else None
        self.independent_judge = IndependentDeterministicJudge()
        self.function_retention_agent = FunctionRetentionAgent(function_retention_backend)
        self.artifact_store = artifact_store

    def run(self, request: WorkflowRequest) -> WorkflowResult:
        from ..agents.quality import WorkflowReassessmentRequired

        self.stage_manifest = []
        try:
            result = self._run(request)
        except WorkflowReassessmentRequired as exc:
            self._complete_diagnostics(request, exc.reasons)
            raise
        except Exception as exc:
            self.stage_manifest.append({"stage": "workflow", "status": "FAILED", "error_type": type(exc).__name__})
            if self.artifact_store:
                run_id = f"{request.target.agi}-{request.herbicide.name}".replace(",", "").replace(" ", "-")
                self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)
            if any(s["stage"] == "input_validation" for s in self.stage_manifest):
                self._complete_diagnostics(request, ["Workflow execution failed: " + type(exc).__name__])
            raise
        self._complete_diagnostics(request, [], result)
        return result

    def _complete_diagnostics(self, request, reasons, result=None):
        from ..agents.learning_agent import LearningAgent
        from ..agents.workflow_review import LEGEND
        from .learning_workflow import LearningWorkflow

        if not (self.require_live_gates or self.diagnostic_review_agent or self.assays):
            return
        run_id = f"{request.target.agi}-{request.herbicide.name}".replace(",", "").replace(" ", "-")

        def record(name, status, detail):
            self.stage_manifest = [s for s in self.stage_manifest if s["stage"] != name]
            self.stage_manifest.append({"stage": name, "status": status, "detail": detail})

        learning = {"status": "PENDING_ASSAYS", "assay_count": len(self.assays), "model_updated": False,
                    "reason": "Real provenance-bearing assays and reviewed candidates are required"}
        if self.assays:
            if self.require_live_gates and any(
                    not assay.provenance or any(p.evidence_type == "synthetic" or p.source.startswith(("mock:", "fixture:"))
                                                for p in assay.provenance) for assay in self.assays):
                learning.update(status="FAILED", reason="Synthetic or unprovenanced assays are not accepted in live workflows")
                reasons = reasons + ["Live assay provenance validation failed"]
            elif result is None:
                learning.update(status="BLOCKED", reason="No reviewed candidates for assay validation")
            elif self.artifact_store is None:
                learning.update(status="BLOCKED", reason="Learning requires a persistent artifact store")
            else:
                try:
                    report, _ = LearningWorkflow(LearningAgent(self.registry), self.artifact_store).run(
                        request, result, self.assays)
                    learning.update(status="COMPLETED", model_updated=True, report=asdict(report))
                except ValueError as exc:
                    learning.update(status="FAILED", error_type=type(exc).__name__)
                    reasons = reasons + ["Assay validation failed"]
        record("learning", learning["status"], learning)
        packet = {"target_agi": request.target.agi, "herbicide": request.herbicide.name,
                  "scope": "workflow_diagnostics_only", "blocking_reasons": list(reasons),
                  "stages": list(self.stage_manifest), "candidate_count": len(result.packets) if result else 0,
                  "candidate_approval": False, "legend": LEGEND}
        record("workflow_review_packet", "COMPLETED", {"scope": "workflow_diagnostics_only"})
        review = self.diagnostic_review_agent.run(packet) if self.diagnostic_review_agent else None
        for role, name in (("review", "workflow_model_review"), ("judge", "workflow_llm_judge")):
            value = review[role] if review else {"status": "UNAVAILABLE", "reason": "No diagnostic model configured"}
            record(name, value["status"], value)
        record("workflow_deterministic_evaluation", "COMPLETED", {
            "approval_permitted": False, "reasons": reasons, "scope": "workflow_diagnostics_only"})
        if self.artifact_store:
            self.artifact_store.write_manifest(run_id, "workflow_review_packet.json", packet)
            self.artifact_store.write_manifest(run_id, "workflow_review.json", review)
            self.artifact_store.write_manifest(run_id, "learning_readiness.json", learning)
            self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)

    def _run(self, request: WorkflowRequest) -> WorkflowResult:
        from ..agents.quality import WorkflowReassessmentRequired

        self.stage_manifest = []
        entry = validate_request(request, self.registry)
        run_id = f"{entry.agi}-{entry.herbicide}".replace(",", "").replace(" ", "-")

        def stage(name, status, detail=None):
            self.stage_manifest.append({"stage": name, "status": status, "detail": detail})
            if self.artifact_store:
                self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)

        def block_remaining(reason):
            recorded = {item["stage"] for item in self.stage_manifest}
            for name in ("pose_ensemble_and_fingerprint", "independent_docking", "pocket_evidence",
                         "constrained_mutation", "candidate_oracles", "scoring_and_review_packets",
                         "model_review", "llm_judge", "deterministic_evaluation"):
                if name not in recorded:
                    stage(name, "NOT_RUN", reason)
            stage("learning", "PENDING_ASSAYS", "No real assay data supplied")

        stage("input_validation", "COMPLETED")
        if self.require_live_gates and (self.evidence_agent is None or self.quality_agent is None):
            stage("preflight", "BLOCKED", "Live workflows require evidence and structure-quality agents")
            raise WorkflowReassessmentRequired(["Missing live evidence/quality gates"])
        blocking_reasons = []
        try:
            evidence = self.evidence_agent.synthesize(entry) if self.evidence_agent else {"facts": []}
            stage("evidence", "COMPLETED" if self.evidence_agent else "UNAVAILABLE", evidence)
        except (RuntimeError, ValueError, OSError) as exc:
            evidence = {"facts": [], "error_type": type(exc).__name__}
            stage("evidence", "UNAVAILABLE", evidence)
            blocking_reasons.append("Evidence synthesis unavailable or invalid; no substitute model used")
        retrieval = getattr(self.evidence_agent, "last_retrieval", None)
        if isinstance(retrieval, dict):
            stage("literature_retrieval", retrieval["status"], retrieval)
        structures = self.structure_backend.predict_ensemble(request.target, entry)
        stage("structure_prediction", "COMPLETED", {"models": len(structures)})
        if not structures or len({item.model_id for item in structures}) != len(structures):
            raise ValueError("A nonempty structure ensemble with unique model IDs is required")
        if any(item.target_agi != request.target.agi for item in structures):
            raise ValueError("Structure ensemble target does not match the request")
        context_errors = [
            error
            for structure in structures
            for error in (
                self.complex_backend.validate_context(structure, entry)
                + validate_structure_context(structure, entry)
            )
        ]
        if self.quality_agent:
            quality = self.quality_agent.assess(request, structures)
            if context_errors:
                quality["status"] = "NEEDS_REASSESSMENT"
                quality["reasons"] = sorted(set(quality["reasons"] + context_errors))
            stage("structure_quality", quality["status"], quality)
            if quality["status"] != "PASSED":
                blocking_reasons.extend(quality["reasons"])
        elif context_errors:
            raise ValueError(f"Structure context validation failed: {context_errors}")
        native_poses = [
            pose
            for ligand in request.native_ligands
            for pose in self.docking_backend.dock(structures, ligand)
        ]
        herbicide_poses = self.docking_backend.dock(structures, request.herbicide)
        model_ids = {item.model_id for item in structures}
        for ligand in [request.herbicide] + request.native_ligands:
            poses = [pose for pose in native_poses + herbicide_poses if pose.ligand_name == ligand.name]
            if not poses or any(pose.model_id not in model_ids for pose in poses):
                raise ValueError(f"Missing or mismatched pose evidence for {ligand.name}")
        stage("pose_ensemble", "COMPLETED", {"native_poses": len(native_poses), "herbicide_poses": len(herbicide_poses),
                                            "scope": "Supplied/cofolded poses; independent docking reported separately"})
        fingerprint = self.fingerprint_agent.compare(
            entry.agi, native_poses, herbicide_poses, request.protected_residues
        )
        stage("pose_ensemble_and_fingerprint", "COMPLETED")
        stage("interaction_fingerprint", "COMPLETED", {"herbicide_selective_positions": fingerprint.herbicide_selective_mutable,
                                                      "shared_protected_positions": fingerprint.shared_protected})
        if self.artifact_store:
            self.artifact_store.write_manifest(run_id, "structure_ensemble.json", structures)
            self.artifact_store.write_manifest(run_id, "pose_ensemble.json", native_poses + herbicide_poses)
            self.artifact_store.write_manifest(run_id, "interaction_fingerprint.json", fingerprint)
        diagnostics = None
        try:
            diagnostics = self.pose_ensemble_agent.run(request, structures) if self.pose_ensemble_agent else None
            stage("independent_docking", "COMPLETED_DIAGNOSTIC_ONLY" if diagnostics else "UNAVAILABLE")
        except (RuntimeError, ValueError, OSError) as exc:
            stage("independent_docking", "FAILED", {"error_type": type(exc).__name__})
            blocking_reasons.append("Independent docking diagnostic failed")
        if self.quality_agent:
            pocket = self.quality_agent.assess_pocket(request, fingerprint, native_poses, herbicide_poses, diagnostics=diagnostics)
            stage("pocket_evidence", pocket["status"], pocket)
            if pocket["status"] != "PASSED":
                blocking_reasons.extend(pocket["reasons"])
        if self.require_live_gates:
            from ..evidence_system.calibration import validate_live_registry
            calibration_reasons = validate_live_registry(
                self.evidence_calibration, target_id=request.target.agi,
                herbicide=request.herbicide.name,
                native_ligands=[ligand.name for ligand in request.native_ligands],
                epsps="epsps" in (str(entry.protein_name) if hasattr(entry, "protein_name") else str(entry)).lower()
                or request.herbicide.name.lower() == "glyphosate",
            )
            stage("calibrated_evidence_baseline", "BLOCKED" if calibration_reasons else "PASSED",
                  {"reasons": calibration_reasons,
                   "calibration_digest": (self.evidence_calibration or {}).get("digest")})
            blocking_reasons.extend(calibration_reasons)
        if self.diagnostics_only:
            blocking_reasons.append("Diagnostic-only run: mutation generation deliberately disabled")
        if blocking_reasons:
            stage("mutation_and_downstream", "BLOCKED", "Diagnostics complete; scientific approval requirements unchanged")
            stage("constrained_mutation", "BLOCKED", {"reasons": blocking_reasons, "candidates": 0})
            stage("candidate_oracles", "BLOCKED", "No approved candidate inputs; mutant-specific artifacts required")
            block_remaining("Candidate design blocked; workflow-level review still runs")
            raise WorkflowReassessmentRequired(sorted(set(blocking_reasons)))
        candidates, rejected = self.mutation_agent.propose(request.target, request.protected_residues, fingerprint)
        stage("constrained_mutation", "COMPLETED", {"candidates": len(candidates)})
        if self.require_live_gates and candidates and self.candidate_evidence_backend is None:
            stage("candidate_oracles", "BLOCKED", "No mutant-specific oracle backend")
            block_remaining("Missing mutant-specific evidence")
            raise WorkflowReassessmentRequired(["Mutant-specific scientific evidence is required"])
        if self.candidate_evidence_backend and candidates:
            try:
                candidates = self.candidate_evidence_backend.evaluate_candidates(request, candidates)
            except WorkflowReassessmentRequired:
                stage("candidate_oracles", "BLOCKED", "Missing protocol-matched mutant evidence")
                block_remaining("Missing mutant-specific evidence")
                raise
            stage("candidate_oracles", "COMPLETED")
        facts = [
            f"Fixed registry pairing: {entry.agi} / {entry.herbicide}.",
            f"Mechanism class: {entry.mechanism_class}.",
            f"Required workflow: {entry.modeling_workflow}.",
        ] + list(evidence.get("facts", []))
        score_packets = [
            self.scoring_agent.score(candidate, native_poses, herbicide_poses, fingerprint)
            for candidate in candidates
        ]
        ranking = pareto_rank(candidates, score_packets)
        function_retention = self.function_retention_agent.assess(
            request.target, candidates, score_packets, request.herbicide, request.native_ligands, structures
        )
        retention_by_mutation = {record.mutation: record for record in function_retention}
        packets = [self.review_agent.review(candidate, scores, facts, [
            "Deterministic function-retention assessment: " + json.dumps(
                asdict(retention_by_mutation[candidate.mutation]), allow_nan=False, sort_keys=True)
        ]) for candidate, scores in zip(candidates, score_packets)]
        stage("scoring_and_review_packets", "COMPLETED")
        stage("model_review", "COMPLETED" if self.review_agent.reasoner else "UNAVAILABLE")
        gated_packets = []
        for packet in packets:
            record = retention_by_mutation[packet.candidate.mutation]
            if record.decision != "MEETS_COMPUTATIONAL_SCREEN":
                failed = record.decision == "FAILS_COMPUTATIONAL_SCREEN"
                packet = replace(
                    packet, recommendation="reject" if failed else "more_computation",
                    status="REJECTED" if failed else "NEEDS_MORE_COMPUTATION",
                    unresolved_uncertainty=packet.unresolved_uncertainty + [
                        f"Function-retention gate: {record.decision}; "
                        f"missing={record.missing_evidence}; failed={record.failed_checks}"
                    ],
                )
            gated_packets.append(packet)
        packets = gated_packets
        failures_by_mutation = {
            packet.candidate.mutation: validate_evaluation_packet(
                packet, request.target.sequence, request.protected_residues, entry
            )
            for packet in packets
        }
        validation_failures = [failure for failures in failures_by_mutation.values() for failure in failures]
        if validation_failures:
            raise ValueError(f"Evaluation packet validation failed: {validation_failures}")
        judge_results = []
        for packet in packets:
            if self.domain_judge:
                judge_results.append(self.domain_judge.judge(packet))
            judge_results.append(self.independent_judge.judge(packet, failures_by_mutation[packet.candidate.mutation]))
        stage("llm_judge", "COMPLETED" if self.domain_judge else "UNAVAILABLE")
        stage("deterministic_evaluation", "COMPLETED")
        stage("learning", "PENDING_ASSAYS", "No learning without provenance-bearing experimental measurements")
        result = WorkflowResult(
            entry, structures, native_poses + herbicide_poses, packets, rejected, [fingerprint], ranking,
            judge_results, function_retention,
        )
        validate_provenance(result)
        if self.artifact_store:
            run_id = f"{entry.agi}-{entry.herbicide}".replace(",", "").replace(" ", "-")
            self.artifact_store.write_manifest(run_id, "structure_ensemble.json", structures)
            self.artifact_store.write_manifest(run_id, "pose_ensemble.json", native_poses + herbicide_poses)
            self.artifact_store.write_manifest(run_id, "interaction_fingerprint.json", fingerprint)
            self.artifact_store.write_manifest(run_id, "evaluation_packets.json", packets)
            self.artifact_store.write_manifest(run_id, "candidate_mutations.json", candidates)
            self.artifact_store.write_manifest(run_id, "rejected_mutations.json", rejected)
            self.artifact_store.write_manifest(run_id, "score_packets.json", score_packets)
            self.artifact_store.write_manifest(run_id, "pareto_ranking.json", ranking)
            self.artifact_store.write_manifest(run_id, "deterministic_validation.json", failures_by_mutation)
            self.artifact_store.write_manifest(run_id, "judge_results.json", judge_results)
            self.artifact_store.write_manifest(run_id, "function_retention_report.json", function_retention)
            self.artifact_store.write_manifest(run_id, "function_retention_thresholds.json", self.function_retention_agent.thresholds)
            self.artifact_store.write_text(run_id, "function_retention_legend.md", RETENTION_LEGEND + "\n")
            self.artifact_store.write_text(
                run_id, "function_retention_table.csv",
                retention_csv(function_retention, request.herbicide.name, [item.name for item in request.native_ligands]),
            )
            self.artifact_store.write_text(
                run_id, "function_retention_table.md",
                retention_markdown(
                    function_retention, request.herbicide.name, [item.name for item in request.native_ligands]
                ),
            )
            renderer = MolstarArtifactRenderer(self.artifact_store)
            visualizations = [
                renderer.render(run_id, structure, index)
                for index, structure in enumerate(structures, start=1)
                if structure.artifact_path
            ]
            self.artifact_store.write_manifest(run_id, "visualizations.json", visualizations)
        return result
