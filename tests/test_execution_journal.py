"""Synthetic software tests for actual execution tracking, not scientific validation."""
import copy
import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock

from herbicide_desensitization_agent.app.orchestrator.execution import ExecutionJournal, payload_hash, validate_journal
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.examples.run_all import make_request


class JournalTests(unittest.TestCase):
    def test_actual_payload_and_file_hashes_persist(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            source, dest = root / "input.dat", root / "output.dat"
            source.write_bytes(b"input bytes")
            journal = ExecutionJournal(root / "execution.json", run_id="test-only")
            def work():
                self.assertEqual(json.loads(journal.path.read_text())["stages"][0]["execution_status"], "RUNNING")
                dest.write_bytes(source.read_bytes().upper())
                return {"file": str(dest), "measurement": None}
            result = journal.execute("structure", work, inputs={"scope": "test"}, input_artifacts=[source], output_artifacts=[dest], evidence_status="INSUFFICIENT_EVIDENCE")
            row = journal.latest("structure")
            self.assertEqual(row["execution_status"], "COMPLETED")
            self.assertEqual(row["evidence_status"], "INSUFFICIENT_EVIDENCE")
            self.assertEqual(row["inputs_sha256"], payload_hash({"scope": "test"}))
            self.assertEqual(row["outputs_sha256"], payload_hash(result))
            self.assertEqual(row["input_artifacts"][0]["sha256"], hashlib.sha256(b"input bytes").hexdigest())
            self.assertEqual(row["output_artifacts"][0]["sha256"], hashlib.sha256(b"INPUT BYTES").hexdigest())
            self.assertIsNotNone(row["started_at"])
            self.assertGreaterEqual(row["runtime_seconds"], 0)
            journal.finalize()
            validate_journal(json.loads(journal.path.read_text()))
            with self.assertRaises(FileExistsError):
                ExecutionJournal(journal.path)

    def test_failures_sanitized_and_independent_stage_runs(self):
        journal = ExecutionJournal()
        def fail():
            raise RuntimeError("Bearer secret-do-not-log")
        with self.assertRaises(RuntimeError):
            journal.execute("docking", fail, inputs={"ligand": "test"})
        self.assertEqual(journal.execute("independent", lambda: 42), 42)
        data = journal.finalize()
        self.assertEqual(data["status"], "COMPLETED_WITH_FAILURES")
        self.assertNotIn("secret-do-not-log", json.dumps(data))
        self.assertEqual(journal.latest("docking")["error_type"], "RuntimeError")
        validate_journal(data)

    def test_skip_and_import_have_no_originating_model_timing(self):
        with tempfile.TemporaryDirectory() as root:
            source = Path(root) / "precomputed.json"
            source.write_text("{}")
            journal = ExecutionJournal()
            for status in sorted(journal.SKIP_STATES):
                journal.skip(status, status, "Explicit test-only reason")
            row = journal.import_artifacts("previous_models", [source])
            self.assertEqual(row["execution_status"], "IMPORTED")
            self.assertIsNone(row["started_at"])
            self.assertIsNone(row["runtime_seconds"])
            validate_journal(journal.finalize())
            with self.assertRaises(ValueError):
                ExecutionJournal().skip("fake", "COMPLETED", "not actually run")

    def test_declared_missing_output_is_failed_not_ready(self):
        journal = ExecutionJournal()
        with self.assertRaises(FileNotFoundError):
            journal.execute("structure", lambda: {}, output_artifacts=["/nonexistent/execution-test-output.pdb"])
        self.assertEqual(journal.latest("structure")["execution_status"], "FAILED")
        validate_journal(journal.finalize())

    def test_repeated_stages_retain_distinct_events(self):
        journal = ExecutionJournal()
        journal.execute("replicate", lambda: 1, inputs={"seed": 1})
        journal.execute("replicate", lambda: 2, inputs={"seed": 2})
        self.assertNotEqual(journal.stages[0]["id"], journal.stages[1]["id"])
        self.assertEqual(len(journal.stages), 2)
        validate_journal(journal.finalize())

    def test_export_validator_rejects_invented_completion(self):
        journal = ExecutionJournal()
        journal.skip("not_executed", "SKIPPED_DEPENDENCY", "missing structures")
        data = copy.deepcopy(journal.finalize())
        data["stages"][0]["execution_status"] = "COMPLETED"
        with self.assertRaises(ValueError):
            validate_journal(data)


class WorkflowExecutionTests(unittest.TestCase):
    def setUp(self):
        self.registry = TargetRegistry.default()
        self.request = make_request(self.registry.entries[0])
        self.backend = MockScientificBackend()

    def workflow(self, **kwargs):
        return WorkflowOrchestrator(self.registry, *([self.backend]*5), **kwargs)

    def test_nomination_with_zero_candidates_does_not_claim_model_work(self):
        workflow = self.workflow()
        workflow.mutation_agent = Mock()
        workflow.mutation_agent.propose.return_value = ([], [{"mutation": "none", "reason": "no eligible contacts"}])
        result = workflow.run(self.request)
        self.assertEqual(result.packets, [])
        for name in ("candidate_oracles", "affinity_scoring", "model_review", "llm_judge", "deterministic_evaluation"):
            row = workflow.execution_journal.latest(name)
            self.assertEqual(row["execution_status"], "SKIPPED_NO_CANDIDATES", name)
            self.assertIsNone(row["runtime_seconds"])
        self.assertEqual(workflow.execution_journal.latest("constrained_mutation")["execution_status"], "COMPLETED")
        self.assertEqual(workflow.execution_journal.latest("function_retention")["execution_status"], "COMPLETED")
        validate_journal(workflow.execution_journal.data)

    def test_ligand_failure_does_not_prevent_other_ligands_or_independent_docking(self):
        original = self.backend.dock
        calls = []
        def dock(structures, ligand):
            calls.append(ligand.name)
            if ligand is self.request.native_ligands[0]:
                raise RuntimeError("native tool failed")
            return original(structures, ligand)
        self.backend.dock = dock
        independent = Mock()
        independent.run.return_value = {"real_call_test_fixture": True}
        workflow = self.workflow(pose_ensemble_agent=independent)
        with self.assertRaises(RuntimeError):
            workflow.run(self.request)
        self.assertEqual(set(calls), {self.request.herbicide.name, *(x.name for x in self.request.native_ligands)})
        independent.run.assert_called_once()
        self.assertEqual(workflow.execution_journal.latest("independent_docking")["execution_status"], "COMPLETED")
        self.assertEqual(workflow.execution_journal.latest("interaction_fingerprint")["execution_status"], "SKIPPED_DEPENDENCY")
        self.assertEqual(workflow.execution_journal.latest("workflow_deterministic_evaluation")["execution_status"], "COMPLETED")
        validate_journal(workflow.execution_journal.data)

    def test_quality_exception_does_not_prevent_independent_docking_or_review(self):
        quality, independent = Mock(), Mock()
        quality.assess.side_effect = RuntimeError("quality crashed")
        quality.assess_pocket.return_value = {"status": "PASSED", "reasons": []}
        independent.run.return_value = {"fixture": True}
        workflow = self.workflow(quality_agent=quality, pose_ensemble_agent=independent)
        from herbicide_desensitization_agent.app.agents.quality import WorkflowReassessmentRequired
        with self.assertRaises(WorkflowReassessmentRequired):
            workflow.run(self.request)
        independent.run.assert_called_once()
        self.assertEqual(workflow.execution_journal.latest("structure_quality")["execution_status"], "FAILED")
        self.assertEqual(workflow.execution_journal.latest("function_retention")["execution_status"], "COMPLETED")
        self.assertEqual(workflow.execution_journal.latest("workflow_deterministic_evaluation")["execution_status"], "COMPLETED")

    def test_model_failure_preserves_deterministic_candidate_assessment(self):
        self.backend.review = Mock(side_effect=RuntimeError("model transport failed"))
        self.backend.judge = Mock(side_effect=RuntimeError("judge transport failed"))
        workflow = self.workflow()
        result = workflow.run(self.request)
        self.assertTrue(result.packets)
        self.assertTrue(result.judge_results)
        self.assertEqual(workflow.execution_journal.latest("model_review")["execution_status"], "FAILED")
        self.assertEqual(workflow.execution_journal.latest("llm_judge")["execution_status"], "FAILED")
        self.assertEqual(workflow.execution_journal.latest("deterministic_evaluation")["execution_status"], "COMPLETED")
        validate_journal(workflow.execution_journal.data)

    def test_external_journal_stays_open_for_fresh_campaign_stages(self):
        journal = ExecutionJournal()
        workflow = self.workflow(journal=journal)
        workflow.run(self.request)
        self.assertEqual(journal.data["status"], "RUNNING")
        journal.execute("dashboard", lambda: {"generated": True})
        validate_journal(journal.finalize())


if __name__ == "__main__":
    unittest.main()
