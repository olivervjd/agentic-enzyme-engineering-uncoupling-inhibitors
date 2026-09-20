from __future__ import annotations

import json
import uuid
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
from .execution import ExecutionJournal


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
        journal: ExecutionJournal | None = None,
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
        self.execution_journal = journal
        self._external_journal = journal is not None
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
        run_id = f"{request.target.agi}-{request.herbicide.name}".replace(",", "").replace(" ", "-")
        if not self._external_journal:
            path = (self.artifact_store.run_directory(run_id) / f"execution-{uuid.uuid4().hex}.json") if self.artifact_store else None
            self.execution_journal = ExecutionJournal(path, run_id=run_id)
        try:
            result = self._run(request)
        except WorkflowReassessmentRequired as exc:
            self._complete_diagnostics(request, exc.reasons)
            raise
        except Exception as exc:
            self._record("workflow", "FAILED", {"error_type": type(exc).__name__})
            if self.artifact_store:
                run_id = f"{request.target.agi}-{request.herbicide.name}".replace(",", "").replace(" ", "-")
                self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)
            if any(s["stage"] == "input_validation" for s in self.stage_manifest):
                self._complete_diagnostics(request, ["Workflow execution failed: " + type(exc).__name__])
            raise
        else:
            self._complete_diagnostics(request, [], result)
            return result
        finally:
            self._record_unvisited()
            if not self._external_journal:
                self.execution_journal.finalize()
            if self.artifact_store:
                self.artifact_store.write_manifest(run_id, "execution_journal.json", self.execution_journal.data)
                self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)

    def _execute(self, name, function, *, inputs=None, dependencies=(), evidence_status="AVAILABLE", **kwargs):
        return self.execution_journal.execute(name, function, inputs=inputs, dependencies=dependencies,
                                              evidence_status=evidence_status, **kwargs)

    def _skip(self, name, status, reason, *, dependencies=(), inputs=None):
        return self.execution_journal.skip(name, status, reason, dependencies=dependencies, inputs=inputs)

    def _record(self, name, status, detail=None):
        event = self.execution_journal.latest(name)
        if event is None:
            skip = "SKIPPED_UNCONFIGURED" if status == "UNAVAILABLE" else "SKIPPED_NO_INPUT" if status == "PENDING_ASSAYS" else "SKIPPED_DEPENDENCY"
            event = self._skip(name, skip, detail if isinstance(detail, str) else "No operation executed for this summary/dependent stage")
        row = {"stage": name, "status": status, "detail": detail,
               **{key: event[key] for key in ("execution_status", "evidence_status", "started_at", "finished_at", "runtime_seconds", "inputs_sha256", "outputs_sha256", "input_artifacts", "output_artifacts", "execution_mode")},
               "execution_event_id": event["id"]}
        self.stage_manifest = [s for s in self.stage_manifest if s["stage"] != name]
        self.stage_manifest.append(row)

    def _record_unvisited(self):
        names = ("input_validation", "evidence", "structure_prediction", "structure_context_validation", "structure_quality",
                 "pose_ensemble", "interaction_fingerprint", "independent_docking", "pocket_evidence",
                 "calibrated_evidence_baseline", "constrained_mutation", "candidate_oracles", "affinity_scoring",
                 "function_retention", "scoring_and_review_packets", "model_review", "llm_judge",
                 "deterministic_evaluation", "learning", "workflow_review_packet", "workflow_model_review",
                 "workflow_llm_judge", "workflow_deterministic_evaluation")
        for name in names:
            if self.execution_journal.latest(name) is None:
                self._skip(name, "SKIPPED_DEPENDENCY", "Required upstream outputs were not produced")
            if not any(r["stage"] == name for r in self.stage_manifest):
                event = self.execution_journal.latest(name)
                self._record(name, event["execution_status"], event.get("reason"))

    def _complete_diagnostics(self, request, reasons, result=None):
        from ..agents.learning_agent import LearningAgent
        from ..agents.workflow_review import LEGEND
        from .learning_workflow import LearningWorkflow

        run_id = f"{request.target.agi}-{request.herbicide.name}".replace(",", "").replace(" ", "-")

        def record(name, status, detail):
            self._record(name, status, detail)

        learning = {"status": "PENDING_ASSAYS", "assay_count": len(self.assays), "model_updated": False,
                    "reason": "Real provenance-bearing assays and reviewed candidates are required"}
        if self.assays:
            if self.require_live_gates and any(
                    not assay.provenance or any(p.evidence_type == "synthetic" or p.source.startswith(("mock:", "fixture:"))
                                                for p in assay.provenance) for assay in self.assays):
                learning.update(status="FAILED", reason="Synthetic or unprovenanced assays are not accepted in live workflows")
                reasons = reasons + ["Live assay provenance validation failed"]
                self._execute("assay_validation", lambda: {"valid": False, "reason": learning["reason"]},
                              inputs=self.assays, evidence_status="FAILED_VALIDATION")
            elif result is None:
                learning.update(status="BLOCKED", reason="No reviewed candidates for assay validation")
            elif self.artifact_store is None:
                learning.update(status="BLOCKED", reason="Learning requires a persistent artifact store")
            else:
                try:
                    # Journal the serializable report; storage side effects occur inside this actual call.
                    def learn():
                        report, _ = LearningWorkflow(LearningAgent(self.registry), self.artifact_store).run(request, result, self.assays)
                        return report
                    report = self._execute("learning", learn, inputs={"request": request, "result": result, "assays": self.assays})
                    learning.update(status="COMPLETED", model_updated=True, report=asdict(report))
                except ValueError as exc:
                    learning.update(status="FAILED", error_type=type(exc).__name__)
                    reasons = reasons + ["Assay validation failed"]
        if self.execution_journal.latest("learning") is None:
            skip = "SKIPPED_NO_INPUT" if not self.assays else "SKIPPED_NO_CANDIDATES" if result is None or not result.packets else "SKIPPED_DEPENDENCY"
            self._skip("learning", skip, learning["reason"], inputs={"assay_count": len(self.assays)})
        record("learning", learning["status"], learning)
        packet = self._execute("workflow_review_packet", lambda: {"target_agi": request.target.agi, "herbicide": request.herbicide.name,
                  "scope": "workflow_diagnostics_only", "blocking_reasons": list(reasons),
                  "stages": list(self.stage_manifest), "candidate_count": len(result.packets) if result else 0,
                  "candidate_approval": False, "legend": LEGEND}, inputs={"reasons": reasons, "stage_events": self.stage_manifest})
        record("workflow_review_packet", "COMPLETED", {"scope": "workflow_diagnostics_only"})
        review = None
        if self.diagnostic_review_agent:
            from ..agents.workflow_review import WorkflowReviewAgent
            try:
                if isinstance(self.diagnostic_review_agent, WorkflowReviewAgent):
                    review = self.diagnostic_review_agent.run(packet, journal=self.execution_journal)
                else:
                    review = self._execute("workflow_review_agents", lambda: self.diagnostic_review_agent.run(packet), inputs=packet,
                                           evidence_status="NOT_ASSESSED")
            except Exception as exc:
                # Independent deterministic review/persistence must survive failed model transport.
                review = {role: {"status": "FAILED", "error_type": type(exc).__name__} for role in ("review", "judge")}
        else:
            for name in ("workflow_model_review", "workflow_llm_judge"):
                self._skip(name, "SKIPPED_UNCONFIGURED", "No diagnostic model configured")
        for role, name in (("review", "workflow_model_review"), ("judge", "workflow_llm_judge")):
            value = review[role] if review else {"status": "UNAVAILABLE", "reason": "No diagnostic model configured"}
            record(name, value["status"], value)
        deterministic = self._execute("workflow_deterministic_evaluation", lambda: {
            "approval_permitted": False, "reasons": list(reasons), "scope": "workflow_diagnostics_only",
            "candidate_approval_claim_absent": packet["candidate_approval"] is False,
            "failed_execution_stages": [e["stage"] for e in self.execution_journal.stages if e["execution_status"] == "FAILED"]},
            inputs={"packet": packet, "review": review}, evidence_status="INSUFFICIENT_EVIDENCE" if reasons else "AVAILABLE")
        record("workflow_deterministic_evaluation", "COMPLETED", deterministic)
        if self.artifact_store:
            self.artifact_store.write_manifest(run_id, "workflow_review_packet.json", packet)
            self.artifact_store.write_manifest(run_id, "workflow_review.json", review)
            self.artifact_store.write_manifest(run_id, "learning_readiness.json", learning)
            self.artifact_store.write_manifest(run_id, "workflow_stages.json", self.stage_manifest)

    def _run(self, request: WorkflowRequest) -> WorkflowResult:
        from ..agents.quality import WorkflowReassessmentRequired

        self.stage_manifest = []
        entry = self._execute("input_validation", lambda: validate_request(request, self.registry), inputs=request, evidence_status="VALIDATED")
        run_id = f"{entry.agi}-{entry.herbicide}".replace(",", "").replace(" ", "-")

        def stage(name, status, detail=None):
            self._record(name, status, detail)
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
        blocking_reasons = []
        if self.require_live_gates and (self.evidence_agent is None or self.quality_agent is None):
            stage("preflight", "BLOCKED", "Live workflows require evidence and structure-quality agents")
            blocking_reasons.append("Missing live evidence/quality gates")
        try:
            evidence = self._execute("evidence", lambda: self.evidence_agent.synthesize(entry), inputs=entry) if self.evidence_agent else {"facts": []}
            stage("evidence", "COMPLETED" if self.evidence_agent else "UNAVAILABLE", evidence)
        except (RuntimeError, ValueError, OSError) as exc:
            evidence = {"facts": [], "error_type": type(exc).__name__}
            stage("evidence", "UNAVAILABLE", evidence)
            blocking_reasons.append("Evidence synthesis unavailable or invalid; no substitute model used")
        retrieval = getattr(self.evidence_agent, "last_retrieval", None)
        if isinstance(retrieval, dict):
            # Retrieval happened within evidence.synthesize; never invent a second runtime.
            stage("literature_retrieval", retrieval["status"], retrieval)
        if self.require_live_gates:
            from ..evidence_system.calibration import validate_live_registry
            calibration_reasons = self._execute("calibrated_evidence_baseline", lambda: validate_live_registry(
                self.evidence_calibration, target_id=request.target.agi, herbicide=request.herbicide.name,
                native_ligands=[ligand.name for ligand in request.native_ligands],
                epsps="epsps" in str(entry).lower() or request.herbicide.name.lower() == "glyphosate"),
                inputs={"calibration": self.evidence_calibration, "request": request},
                evidence_status=lambda reasons: "INSUFFICIENT_EVIDENCE" if reasons else "VALIDATED")
            stage("calibrated_evidence_baseline", "BLOCKED" if calibration_reasons else "PASSED", {"reasons": calibration_reasons})
            blocking_reasons.extend(calibration_reasons)
        else:
            self._skip("calibrated_evidence_baseline", "SKIPPED_UNCONFIGURED", "Synthetic/demo mode does not execute prospective calibration gates")

        def predict():
            values = self.structure_backend.predict_ensemble(request.target, entry)
            if not values or len({item.model_id for item in values}) != len(values):
                raise ValueError("A nonempty structure ensemble with unique model IDs is required")
            if any(item.target_agi != request.target.agi for item in values):
                raise ValueError("Structure ensemble target does not match the request")
            return values
        structures = self._execute("structure_prediction", predict, inputs={"target": request.target, "registry": entry})
        stage("structure_prediction", "COMPLETED", {"models": len(structures)})
        context_errors = self._execute("structure_context_validation", lambda: [
            error
            for structure in structures
            for error in (
                self.complex_backend.validate_context(structure, entry)
                + validate_structure_context(structure, entry)
            )
        ], inputs={"structures": structures, "registry": entry}, evidence_status=lambda errors: "FAILED_VALIDATION" if errors else "VALIDATED")
        deferred_error = None
        if self.quality_agent:
            try:
                quality = self._execute("structure_quality", lambda: self.quality_agent.assess(request, structures),
                                        inputs={"request": request, "structures": structures},
                                        evidence_status=lambda result: "VALIDATED" if result["status"] == "PASSED" and not context_errors else "INSUFFICIENT_EVIDENCE")
                if context_errors:
                    quality["status"] = "NEEDS_REASSESSMENT"
                    quality["reasons"] = sorted(set(quality["reasons"] + context_errors))
                stage("structure_quality", quality["status"], quality)
                if quality["status"] != "PASSED":
                    blocking_reasons.extend(quality["reasons"])
            except (RuntimeError, ValueError, OSError) as exc:
                stage("structure_quality", "FAILED", {"error_type": type(exc).__name__})
                blocking_reasons.append("Structure quality execution failed")
        elif context_errors:
            deferred_error = ValueError(f"Structure context validation failed: {context_errors}")
        else:
            self._skip("structure_quality", "SKIPPED_UNCONFIGURED", "No structure quality agent configured")
        native_poses, herbicide_poses = [], []
        model_ids = {item.model_id for item in structures}
        pose_failures = []
        def assemble_poses():
            for ligand in request.native_ligands + [request.herbicide]:
                try:
                    def dock():
                        poses = self.docking_backend.dock(structures, ligand)
                        if not poses or any(p.model_id not in model_ids or p.ligand_name != ligand.name for p in poses):
                            raise ValueError(f"Missing or mismatched pose evidence for {ligand.name}")
                        return poses
                    poses = self._execute("ligand_placement:" + ligand.name, dock,
                                          inputs={"structures": structures, "ligand": ligand}, dependencies=["structure_prediction"])
                    (herbicide_poses if ligand is request.herbicide else native_poses).extend(poses)
                except (RuntimeError, ValueError, OSError) as exc:
                    pose_failures.append(exc)
            if pose_failures:
                raise pose_failures[0]
            return native_poses + herbicide_poses
        try:
            self._execute("pose_ensemble", assemble_poses, inputs={"structures": structures, "ligands": request.native_ligands + [request.herbicide]})
            stage("pose_ensemble", "COMPLETED", {"native_poses": len(native_poses), "herbicide_poses": len(herbicide_poses),
                                                "scope": "Supplied/cofolded poses; independent docking reported separately"})
        except (RuntimeError, ValueError, OSError) as exc:
            deferred_error = deferred_error or exc
            stage("pose_ensemble", "FAILED", {"error_type": type(exc).__name__})
        diagnostics = None
        try:
            diagnostics = self._execute("independent_docking", lambda: self.pose_ensemble_agent.run(request, structures),
                                        inputs={"request": request, "structures": structures}) if self.pose_ensemble_agent else None
            stage("independent_docking", "COMPLETED_DIAGNOSTIC_ONLY" if diagnostics else "UNAVAILABLE")
        except (RuntimeError, ValueError, OSError) as exc:
            stage("independent_docking", "FAILED", {"error_type": type(exc).__name__})
            blocking_reasons.append("Independent docking diagnostic failed")
        if deferred_error:
            # Literature/calibration, other ligands and independent docking have
            # already run. Contact selection cannot consume partial ligand data.
            raise deferred_error
        fingerprint = self._execute("interaction_fingerprint", lambda: self.fingerprint_agent.compare(
            entry.agi, native_poses, herbicide_poses, request.protected_residues),
            inputs={"native_poses": native_poses, "herbicide_poses": herbicide_poses, "protected": request.protected_residues})
        stage("interaction_fingerprint", "COMPLETED", {"herbicide_selective_positions": fingerprint.herbicide_selective_mutable,
                                                      "shared_protected_positions": fingerprint.shared_protected})
        if self.artifact_store:
            self.artifact_store.write_manifest(run_id, "structure_ensemble.json", structures)
            self.artifact_store.write_manifest(run_id, "pose_ensemble.json", native_poses + herbicide_poses)
            self.artifact_store.write_manifest(run_id, "interaction_fingerprint.json", fingerprint)
        if self.quality_agent:
            try:
                pocket = self._execute("pocket_evidence", lambda: self.quality_agent.assess_pocket(request, fingerprint, native_poses, herbicide_poses, diagnostics=diagnostics),
                                       inputs={"request": request, "fingerprint": fingerprint, "native": native_poses, "herbicide": herbicide_poses, "diagnostics": diagnostics},
                                       evidence_status=lambda result: "VALIDATED" if result["status"] == "PASSED" else "INSUFFICIENT_EVIDENCE")
                stage("pocket_evidence", pocket["status"], pocket)
                if pocket["status"] != "PASSED":
                    blocking_reasons.extend(pocket["reasons"])
            except (RuntimeError, ValueError, OSError) as exc:
                stage("pocket_evidence", "FAILED", {"error_type": type(exc).__name__})
                blocking_reasons.append("Pocket evidence execution failed")
        else:
            self._skip("pocket_evidence", "SKIPPED_UNCONFIGURED", "No pocket quality agent configured")
        if self.diagnostics_only:
            blocking_reasons.append("Diagnostic-only run: mutation generation deliberately disabled")
        if blocking_reasons:
            # WT evidence assessment is useful even when mutation design is
            # scientifically blocked. Do not label candidate-only work complete.
            try:
                wt_retention = self._execute("function_retention", lambda: self.function_retention_agent.assess(
                    request.target, [], [], request.herbicide, request.native_ligands, structures),
                    inputs={"scope": "WT_only_no_mutation_approved", "request": request, "structures": structures},
                    evidence_status=lambda rows: "INSUFFICIENT_EVIDENCE" if any(r.missing_evidence for r in rows) else "AVAILABLE")
                stage("function_retention", "COMPLETED", {"scope": "WT_only", "records": len(wt_retention)})
                if self.artifact_store:
                    self.artifact_store.write_manifest(run_id, "wt_function_retention_report.json", wt_retention)
            except (RuntimeError, ValueError, OSError) as exc:
                stage("function_retention", "FAILED", {"scope": "WT_only", "error_type": type(exc).__name__})
                blocking_reasons.append("WT function-retention evidence execution failed")
            stage("mutation_and_downstream", "BLOCKED", "Diagnostics complete; scientific approval requirements unchanged")
            stage("constrained_mutation", "BLOCKED", {"reasons": blocking_reasons, "candidates": 0})
            stage("candidate_oracles", "BLOCKED", "No approved candidate inputs; mutant-specific artifacts required")
            block_remaining("Candidate design blocked; workflow-level review still runs")
            raise WorkflowReassessmentRequired(sorted(set(blocking_reasons)))
        candidates, rejected = self._execute("constrained_mutation", lambda: self.mutation_agent.propose(request.target, request.protected_residues, fingerprint),
                                             inputs={"target": request.target, "protected": request.protected_residues, "fingerprint": fingerprint})
        stage("constrained_mutation", "COMPLETED", {"candidates": len(candidates)})
        if self.require_live_gates and candidates and self.candidate_evidence_backend is None:
            stage("candidate_oracles", "BLOCKED", "No mutant-specific oracle backend")
            block_remaining("Missing mutant-specific evidence")
            raise WorkflowReassessmentRequired(["Mutant-specific scientific evidence is required"])
        if self.candidate_evidence_backend and candidates:
            try:
                candidates = self._execute("candidate_oracles", lambda: self.candidate_evidence_backend.evaluate_candidates(request, candidates),
                                           inputs={"request": request, "candidates": candidates})
            except WorkflowReassessmentRequired:
                stage("candidate_oracles", "BLOCKED", "Missing protocol-matched mutant evidence")
                block_remaining("Missing mutant-specific evidence")
                raise
            stage("candidate_oracles", "COMPLETED")
        elif not candidates:
            self._skip("candidate_oracles", "SKIPPED_NO_CANDIDATES", "Mutation selection completed with zero eligible candidates")
            stage("candidate_oracles", "NOT_RUN", "No eligible candidates")
        else:
            self._skip("candidate_oracles", "SKIPPED_UNCONFIGURED", "No mutant oracle backend configured in synthetic/demo mode")
            stage("candidate_oracles", "UNAVAILABLE", "No mutant oracle backend")
        facts = [
            f"Fixed registry pairing: {entry.agi} / {entry.herbicide}.",
            f"Mechanism class: {entry.mechanism_class}.",
            f"Required workflow: {entry.modeling_workflow}.",
        ] + list(evidence.get("facts", []))
        if candidates:
            score_packets = self._execute("affinity_scoring", lambda: [
                self.scoring_agent.score(candidate, native_poses, herbicide_poses, fingerprint) for candidate in candidates],
                inputs={"candidates": candidates, "native_poses": native_poses, "herbicide_poses": herbicide_poses, "fingerprint": fingerprint})
        else:
            score_packets = []
            self._skip("affinity_scoring", "SKIPPED_NO_CANDIDATES", "No candidate-specific affinity scoring can run")
        ranking = self._execute("pareto_ranking", lambda: pareto_rank(candidates, score_packets), inputs={"candidates": candidates, "scores": score_packets}) if candidates else []
        function_retention = self._execute("function_retention", lambda: self.function_retention_agent.assess(
            request.target, candidates, score_packets, request.herbicide, request.native_ligands, structures),
            inputs={"target": request.target, "candidates": candidates, "scores": score_packets, "herbicide": request.herbicide, "native": request.native_ligands, "structures": structures},
            evidence_status=lambda records: "INSUFFICIENT_EVIDENCE" if any(r.missing_evidence for r in records) else "AVAILABLE")
        retention_by_mutation = {record.mutation: record for record in function_retention}
        def build_packets():
            values = []
            for candidate, scores in zip(candidates, score_packets):
                predictions = ["Deterministic function-retention assessment: " + json.dumps(
                    asdict(retention_by_mutation[candidate.mutation]), allow_nan=False, sort_keys=True)]
                if self.review_agent.reasoner:
                    try:
                        packet = self._execute("candidate_model_review:" + candidate.mutation,
                            lambda: self.review_agent.review(candidate, scores, facts, predictions),
                            inputs={"candidate": candidate, "scores": scores, "facts": facts, "predictions": predictions}, evidence_status="NOT_ASSESSED")
                    except (RuntimeError, ValueError, OSError):
                        # Preserve scientific packet and deterministic review after model failure.
                        packet = CandidateReviewAgent(None).review(candidate, scores, facts, predictions)
                        packet = replace(packet, unresolved_uncertainty=packet.unresolved_uncertainty + ["Configured model review failed; deterministic evidence packet retained"])
                else:
                    packet = self.review_agent.review(candidate, scores, facts, predictions)
                values.append(packet)
            return values
        if candidates:
            model_event = self.execution_journal.start("model_review", inputs={"candidates": candidates, "scores": score_packets, "facts": facts}) if self.review_agent.reasoner else None
            try:
                packets = self._execute("scoring_and_review_packets", build_packets, inputs={"candidates": candidates, "scores": score_packets, "retention": function_retention})
                if model_event:
                    failed = any(e["execution_status"] == "FAILED" for e in self.execution_journal.stages if e["stage"].startswith("candidate_model_review:"))
                    if failed:
                        self.execution_journal.fail(model_event, RuntimeError("Candidate model review failed"))
                    else:
                        self.execution_journal.finish(model_event, outputs=packets, evidence_status="NOT_ASSESSED")
            except Exception as exc:
                if model_event:
                    self.execution_journal.fail(model_event, exc)
                raise
            stage("scoring_and_review_packets", "COMPLETED")
            if self.review_agent.reasoner:
                # Per-candidate journal events carry genuine model timings.
                stage("model_review", model_event["execution_status"], {"events": [e["id"] for e in self.execution_journal.stages if e["stage"].startswith("candidate_model_review:")]})
            else:
                self._skip("model_review", "SKIPPED_UNCONFIGURED", "No candidate review model configured")
                stage("model_review", "UNAVAILABLE")
        else:
            packets = []
            for name in ("scoring_and_review_packets", "model_review", "llm_judge", "deterministic_evaluation"):
                self._skip(name, "SKIPPED_NO_CANDIDATES", "No selected mutations require candidate review")
                stage(name, "NOT_RUN", "Mutation nomination produced no eligible candidates")
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
        failures_by_mutation = self._execute("candidate_packet_validation", lambda: {
            packet.candidate.mutation: validate_evaluation_packet(
                packet, request.target.sequence, request.protected_residues, entry
            )
            for packet in packets
        }, inputs={"packets": packets, "request": request, "registry": entry}, evidence_status=lambda results: "FAILED_VALIDATION" if any(results.values()) else "VALIDATED") if packets else {}
        validation_failures = [failure for failures in failures_by_mutation.values() for failure in failures]
        if validation_failures:
            raise ValueError(f"Evaluation packet validation failed: {validation_failures}")
        judge_results = []
        model_judge_event = self.execution_journal.start("llm_judge", inputs=packets) if packets and self.domain_judge else None
        deterministic_event = self.execution_journal.start("deterministic_evaluation", inputs={"packets": packets, "failures": failures_by_mutation}) if packets else None
        try:
            for packet in packets:
                if self.domain_judge:
                    try:
                        judge_results.append(self._execute("candidate_llm_judge:" + packet.candidate.mutation, lambda: self.domain_judge.judge(packet), inputs=packet, evidence_status="NOT_ASSESSED"))
                    except (RuntimeError, ValueError, OSError):
                        pass  # Logged failed execution; deterministic judge still runs.
                judge_results.append(self._execute("candidate_deterministic_evaluation:" + packet.candidate.mutation,
                                     lambda: self.independent_judge.judge(packet, failures_by_mutation[packet.candidate.mutation]),
                                     inputs={"packet": packet, "failures": failures_by_mutation[packet.candidate.mutation]}))
            if model_judge_event:
                failed = any(e["execution_status"] == "FAILED" for e in self.execution_journal.stages if e["stage"].startswith("candidate_llm_judge:"))
                if failed:
                    self.execution_journal.fail(model_judge_event, RuntimeError("Candidate model judge failed"))
                else:
                    self.execution_journal.finish(model_judge_event, outputs=[r for r in judge_results if r.judge_id != self.independent_judge.judge_id], evidence_status="NOT_ASSESSED")
            if deterministic_event:
                self.execution_journal.finish(deterministic_event, outputs=[r for r in judge_results if r.judge_id == self.independent_judge.judge_id])
        except Exception as exc:
            for event in (model_judge_event, deterministic_event):
                if event:
                    self.execution_journal.fail(event, exc)
            raise
        if packets:
            stage("llm_judge", model_judge_event["execution_status"] if self.domain_judge else "UNAVAILABLE")
            stage("deterministic_evaluation", "COMPLETED")
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
