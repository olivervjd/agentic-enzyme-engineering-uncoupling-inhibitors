import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from herbicide_desensitization_agent.app.agents.evidence import CitedEvidenceAgent
from herbicide_desensitization_agent.app.agents.quality import WorkflowReassessmentRequired
from herbicide_desensitization_agent.app.agents.workflow_review import WorkflowReviewAgent
from herbicide_desensitization_agent.app.backends.calibrated_epsps import CalibratedEPSPSBackend
from herbicide_desensitization_agent.app.backends.literature import EuropePMCRetriever
from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.backends.openai_json import UnavailableModelTransport, response_schema
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import AssayResult, Provenance
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.registry = TargetRegistry.default()
        self.entry = self.registry.entries[0]
        self.request = make_request(self.entry)
        self.backend = MockScientificBackend()
        self.quality = Mock()
        self.quality.assess.return_value = {"status": "NEEDS_REASSESSMENT", "reasons": ["chemistry"]}
        self.quality.assess_pocket.return_value = {"status": "PASSED", "reasons": []}
        self.evidence = Mock()
        self.evidence.synthesize.return_value = {"facts": []}
        self.review = Mock(model="review-fixture", return_value={
            "summary": "Diagnostic fixture", "issues": ["chemistry"],
            "next_actions": ["validate chemistry"], "recommendation": "revise"})
        self.judge = Mock(model="judge-fixture", return_value=self.review.return_value)

    def workflow(self, **options):
        return WorkflowOrchestrator(self.registry, *([self.backend] * 5),
            artifact_store=ArtifactStore(self.root), quality_agent=self.quality,
            evidence_agent=self.evidence, require_live_gates=True,
            diagnostic_review_agent=WorkflowReviewAgent(self.review, self.judge), **options)

    def test_quality_failure_runs_reviews_without_mutating(self):
        workflow = self.workflow()
        workflow.mutation_agent = Mock()
        with self.assertRaises(WorkflowReassessmentRequired):
            workflow.run(self.request)
        workflow.mutation_agent.propose.assert_not_called()
        self.review.assert_called_once()
        self.judge.assert_called_once()
        self.assertEqual(self.judge.call_args.args[0], "workflow_judge")
        self.assertIn("review_to_audit", self.judge.call_args.args[1])
        self.assertTrue(list(self.root.rglob("interaction_fingerprint.json")))
        self.assertTrue(list(self.root.rglob("workflow_review_packet.json")))
        self.assertEqual(next(s["status"] for s in workflow.stage_manifest if s["stage"] == "learning"), "PENDING_ASSAYS")

    def test_evidence_unavailable_does_not_stop_diagnostics_or_allow_mutation(self):
        self.evidence.synthesize.side_effect = RuntimeError("not accessible")
        self.quality.assess.return_value = {"status": "PASSED", "reasons": []}
        workflow = self.workflow()
        workflow.mutation_agent = Mock()
        with self.assertRaises(WorkflowReassessmentRequired):
            workflow.run(self.request)
        workflow.mutation_agent.propose.assert_not_called()
        self.review.assert_called_once()
        self.assertEqual(next(s["status"] for s in workflow.stage_manifest if s["stage"] == "evidence"), "UNAVAILABLE")

    def test_diagnostic_flag_never_designs_even_with_passed_gates(self):
        self.quality.assess.return_value = {"status": "PASSED", "reasons": []}
        workflow = self.workflow(diagnostics_only=True)
        workflow.mutation_agent = Mock()
        with self.assertRaisesRegex(WorkflowReassessmentRequired, "Diagnostic-only"):
            workflow.run(self.request)
        workflow.mutation_agent.propose.assert_not_called()

    def test_independent_docking_failure_still_reaches_review(self):
        docking = Mock()
        docking.run.side_effect = RuntimeError("fixture failure")
        workflow = self.workflow(pose_ensemble_agent=docking)
        with self.assertRaises(WorkflowReassessmentRequired):
            workflow.run(self.request)
        self.review.assert_called_once()
        self.assertEqual(next(s["status"] for s in workflow.stage_manifest if s["stage"] == "independent_docking"), "FAILED")

    def test_failed_review_does_not_stop_judge_or_leak_errors(self):
        self.review.side_effect = RuntimeError("secret-token-value")
        report = WorkflowReviewAgent(self.review, self.judge).run({"blocking_reasons": ["fixture"]})
        self.assertEqual(report["review"]["status"], "FAILED")
        self.assertEqual(report["judge"]["status"], "COMPLETED")
        self.assertNotIn("secret-token-value", json.dumps(report))
        self.assertFalse(report["candidate_approval"])

    def test_bad_model_output_is_not_reported_completed(self):
        self.review.return_value = {"summary": "fiction", "issues": [], "next_actions": [], "recommendation": "advance"}
        report = WorkflowReviewAgent(self.review).run({"blocking_reasons": []})
        self.assertEqual(report["review"]["status"], "FAILED")

    def test_new_operation_schemas_are_strict(self):
        for name in ("workflow_review", "workflow_judge"):
            schema = response_schema(name)
            self.assertFalse(schema["additionalProperties"])
            self.assertEqual(set(schema["properties"]), set(schema["required"]))

    def test_cached_only_does_not_launch_prediction(self):
        runner = Mock(root=self.root)
        backend = CalibratedEPSPSBackend(runner, cached_only=True)
        entry = self.registry.get("AT2G45300", "glyphosate")
        with self.assertRaisesRegex(ValueError, "no GPU"):
            backend.predict_ensemble(make_request(entry).target, entry)
        runner.run.assert_not_called()

    def test_learning_is_wired_to_reviewed_candidates(self):
        fixture = WorkflowOrchestrator(self.registry, *([self.backend] * 5)).run(self.request)
        assay = AssayResult("fixture-assay", self.entry.agi, self.entry.herbicide,
            fixture.packets[0].candidate.mutation,
            {"herbicide_response_fold_change": 1.0, "protein_expression_fraction": 1.0, "native_activity_fraction": 1.0},
            {"herbicide_response_fold_change": "fold", "protein_expression_fraction": "fraction", "native_activity_fraction": "fraction"},
            [Provenance("fixture://assay", "test-only", "synthetic")])
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 5),
                                        artifact_store=ArtifactStore(self.root), assays=[assay])
        workflow.run(self.request)
        self.assertTrue(list(self.root.rglob("model_recalibration_report.json")))
        self.assertEqual(next(s["status"] for s in workflow.stage_manifest if s["stage"] == "learning"), "COMPLETED")


class LiteratureTests(unittest.TestCase):
    def test_retrieves_deduplicated_abstracts_with_provenance(self):
        record = {"source": "MED", "id": "123", "title": "EPSPS", "abstractText": "Untrusted abstract"}
        raw = json.dumps({"resultList": {"result": [record, record, {"source": "MED", "id": "456"}]}}).encode()
        with tempfile.TemporaryDirectory() as tmp:
            opener = Mock(return_value=io.BytesIO(raw))
            entry = TargetRegistry.default().get("AT2G45300", "glyphosate")
            report = EuropePMCRetriever(tmp, opener=opener).retrieve(entry)
            self.assertEqual(len(report["sources"]), 1)
            self.assertEqual(report["sources"][0]["evidence_level"], "abstract_only_unreviewed")
            self.assertTrue(Path(report["sources"][0]["artifact"]).is_file())
            self.assertIn("resultType=core", opener.call_args.args[0].full_url)
            self.assertNotIn("Authorization", opener.call_args.args[0].headers)

    def test_retrieval_survives_unavailable_rosalind(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "source.txt").write_text("fixture")
            entry = TargetRegistry.default().entries[0]
            (root / "evidence.json").write_text(json.dumps({"target_agi": entry.agi, "sources": [{
                "id": "s", "url": "https://example.org/s", "scope": "fixture", "artifact": "source.txt",
                "sha256": hashlib.sha256(b"fixture").hexdigest()}], "claims": []}))
            retriever = Mock()
            retriever.retrieve.return_value = {"status": "NO_RESULTS", "sources": []}
            agent = CitedEvidenceAgent(root / "evidence.json", UnavailableModelTransport("rosalind", "unavailable"), retriever)
            with self.assertRaises(RuntimeError):
                agent.synthesize(entry)
            self.assertEqual(agent.last_retrieval["status"], "NO_RESULTS")

    def test_invalid_search_limit_rejected(self):
        with self.assertRaises(ValueError):
            EuropePMCRetriever("unused", limit=100)
