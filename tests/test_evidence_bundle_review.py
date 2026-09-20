import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from herbicide_desensitization_agent.examples.review_evidence_bundle import compact_packet, main, run_review


class EvidenceBundleReviewTests(unittest.TestCase):
    def setUp(self):
        self.bundle = {"objective": "Retain native function", "conclusion": {"nominations": "No nominations"},
                       "calibration": {"status": "INSUFFICIENT_EVIDENCE", "thresholds": {"tm": 0.9}},
                       "limitations": ["No biochemical validation"],
                       "decisions": [{"mutation": "A1V", "decision": "INSUFFICIENT_EVIDENCE"}],
                       "tables": {"model_calibration": [{"control": "WT", "pass_fail": "NOT_EVALUATED"}],
                                  "residue_interactions": [{"large": "excluded"}]},
                       "irrelevant_metadata": {"private": "never-transmitted"}}
        self.output = {"summary": "Calibration incomplete", "issues": ["Missing assays"],
                       "next_actions": ["Gather controls"], "recommendation": "revise"}

    def execute(self, factory, key="test-secret"):
        return run_review(self.bundle, api_key=key, credential_source="codex_api_key",
                          source_provenance={"artifact": "evidence_system.json", "sha256": "fixture"},
                          transport_factory=factory)

    def test_calls_are_separate_and_cannot_mutate_scientific_input(self):
        original = copy.deepcopy(self.bundle)
        calls = []
        def factory(model, api_key):
            self.assertEqual(model, "gpt-5.6-luna")
            def transport(operation, payload):
                calls.append((operation, copy.deepcopy(payload)))
                payload["packet"]["calibration"]["thresholds"]["tm"] = 0
                return copy.deepcopy(self.output)
            return transport
        report = self.execute(factory)
        self.assertEqual(self.bundle, original)
        self.assertEqual([item[0] for item in calls], ["workflow_review", "workflow_judge"])
        self.assertEqual(calls[1][1]["packet"]["calibration"]["thresholds"]["tm"], 0.9)
        self.assertEqual(calls[1][1]["review_to_audit"]["status"], "COMPLETED")
        self.assertEqual(report["decision_counts"], {"INSUFFICIENT_EVIDENCE": 1})
        self.assertTrue(report["scientific_conclusions_unchanged"])
        self.assertFalse(report["candidate_approval"])
        self.assertEqual(report["status"], "COMPLETED")
        self.assertNotIn("never-transmitted", json.dumps(calls))
        self.assertNotIn("test-secret", json.dumps(report))

    def test_failure_never_persists_exception_body_or_fabricates_output(self):
        def factory(*args, **kwargs):
            def transport(*args):
                raise RuntimeError("provider body test-secret")
            return transport
        report = self.execute(factory)
        self.assertEqual(report["status"], "INCOMPLETE")
        for role in ("review", "judge"):
            self.assertEqual(report[role]["status"], "FAILED")
            self.assertIsNone(report[role]["output"])
        self.assertNotIn("provider body", json.dumps(report))
        self.assertNotIn("test-secret", json.dumps(report))

    def test_model_cannot_add_decision_or_numerical_output_fields(self):
        output = {**self.output, "decision": "MEETS_COMPUTATIONAL_SCREEN"}
        report = self.execute(lambda *a, **kw: lambda *args: output)
        self.assertEqual(report["status"], "INCOMPLETE")
        self.assertEqual(report["review"]["error_type"], "ValueError")
        self.assertEqual(self.bundle["decisions"][0]["decision"], "INSUFFICIENT_EVIDENCE")

    def test_even_accidentally_echoed_key_is_redacted(self):
        output = {**self.output, "summary": "accidental test-secret echo"}
        report = self.execute(lambda *a, **kw: lambda *args: output)
        self.assertNotIn("test-secret", json.dumps(report))
        self.assertIn("[REDACTED]", report["review"]["output"]["summary"])

    def test_missing_key_does_not_call_transport(self):
        def forbidden(*args, **kwargs):
            self.fail("Transport must not run without credentials")
        report = self.execute(forbidden, key=None)
        self.assertEqual(report["review"]["status"], "UNAVAILABLE")
        self.assertEqual(report["status"], "INCOMPLETE")

    def test_compact_control_rows_are_bounded_and_omission_is_recorded(self):
        self.bundle["tables"]["model_calibration"] *= 60
        packet = compact_packet(self.bundle)
        self.assertEqual(len(packet["model_calibration"]), 50)
        self.assertEqual(packet["control_rows_total"], 60)

    def test_compact_new_control_payload_excludes_atom_coordinates_and_sequences(self):
        self.bundle['tables']['model_calibration']=[{'control':'Vina','predicted_outcome':{
            'top_pose_score_mean':-4.1,'returned_pose_count':15,'independent_docking_seed_count':3,
            'preparation':{'protein_sequence':'SECRET_SEQUENCE_MARKER','protein_heavy_atoms_added_by_templates':[{'xyz':[1,2,3]}]},
            'protein_transform':{'rotation':[[1,0,0]]},'pose_clusters':[{'cluster':0}]}}]
        self.bundle['calibration']['core_gate_results']={'core_calibration':{'specification':{'baseline_evidence':{
            'documented_msa':{'validated':True,'acceptance_calibrated':False,'provenance':{'sequence':'SECRET_SEQUENCE_MARKER'}}}},'wt_replicates':[{'coordinates':'COORDINATE_MARKER'}]}}
        packet=compact_packet(self.bundle); encoded=json.dumps(packet)
        self.assertNotIn('SECRET_SEQUENCE_MARKER',encoded); self.assertNotIn('COORDINATE_MARKER',encoded)
        self.assertNotIn('protein_transform',encoded); self.assertNotIn('xyz',encoded)
        summary=packet['model_calibration'][0]['predicted_outcome']
        self.assertEqual(summary['top_pose_score_mean'],-4.1); self.assertEqual(summary['pose_cluster_count'],1)
        self.assertTrue(packet['calibration']['baseline_evidence']['documented_msa']['validated'])

    def test_cli_preserves_exact_source_bytes_and_writes_separate_failure_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "evidence_system.json"
            original = json.dumps(self.bundle).encode()
            source.write_bytes(original)
            with patch("sys.argv", ["review_evidence_bundle", "--input", str(source)]), \
                    patch("herbicide_desensitization_agent.examples.review_evidence_bundle.resolve_api_key", return_value=(None, "unavailable")):
                self.assertEqual(main(), 2)
            self.assertEqual(source.read_bytes(), original)
            report = json.loads((source.parent / "evidence_model_review.json").read_text())
            self.assertEqual(report["status"], "INCOMPLETE")
            raw_digest = hashlib.sha256(original).hexdigest()
            self.assertEqual(report["source_bundle_sha256"], raw_digest)
            for role in ("review", "judge"):
                self.assertEqual(report[role]["provenance"]["source"]["sha256"], raw_digest)

    def test_external_source_change_during_transport_marks_snapshot_stale(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "evidence_system.json"
            original = json.dumps(self.bundle, indent=3).encode()
            source.write_bytes(original)
            replacement = {**self.bundle, "limitations": ["External updated evidence"]}
            def factory(*args, **kwargs):
                def transport(*args):
                    source.write_text(json.dumps(replacement))
                    return copy.deepcopy(self.output)
                return transport
            def injected_review(bundle, **kwargs):
                return run_review(bundle, **kwargs, transport_factory=factory)
            prefix = "herbicide_desensitization_agent.examples.review_evidence_bundle."
            with patch("sys.argv", ["review_evidence_bundle", "--input", str(source)]), \
                    patch(prefix + "resolve_api_key", return_value=("test-secret", "codex_api_key")), \
                    patch(prefix + "run_review", side_effect=injected_review):
                self.assertEqual(main(), 2)
            report = json.loads((source.parent / "evidence_model_review.json").read_text())
            self.assertEqual(report["status"], "STALE")
            self.assertFalse(report["source_snapshot_current"])
            self.assertEqual(json.loads(source.read_bytes()), replacement)
            old_hash = hashlib.sha256(original).hexdigest()
            self.assertEqual(report["source_bundle_sha256"], old_hash)
            for role in ("review", "judge"):
                self.assertEqual(report[role]["status"], "STALE")
                self.assertEqual(report[role]["execution_status"], "COMPLETED")
                self.assertEqual(report[role]["provenance"]["source"]["sha256"], old_hash)
            self.assertNotIn("test-secret", json.dumps(report))


if __name__ == "__main__":
    unittest.main()
