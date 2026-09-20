import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock, patch
from urllib.error import HTTPError

from herbicide_desensitization_agent.app.agents.evidence import CitedEvidenceAgent
from herbicide_desensitization_agent.app.agents.pipeline_agents import CandidateReviewAgent, MultiOracleScoringAgent
from herbicide_desensitization_agent.app.agents.quality import StructureQualityAgent, WorkflowReassessmentRequired
from herbicide_desensitization_agent.app.agents.pose_ensemble import PoseEnsembleAgent
from herbicide_desensitization_agent.app.backends.calibration import CONTEXTS, SEEDS, validate_a3m, replace_msa_query, validate_record_coverage, write
from herbicide_desensitization_agent.app.backends.calibrated_epsps import ArtifactCandidateOracles
from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.backends.openai_json import OpenAIJSONTransport, response_schema
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import Provenance
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


class LiveIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.registry = TargetRegistry.default()
        self.entry = self.registry.entries[0]
        self.request = make_request(self.entry)
        self.backend = MockScientificBackend()
        self.workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 5))
        self.fixture = self.workflow.run(self.request)

    def test_msa_rejects_query_only(self):
        with self.assertRaisesRegex(ValueError, "query-only"):
            validate_a3m(">q\nACDE\n", "ACDE")

    def test_msa_rejects_wrong_query_and_width(self):
        for text in (">q\nACDF\n>s\nACDE\n", ">q\nACDE\n>s\nACD\n"):
            with self.assertRaises(ValueError):
                validate_a3m(text, "ACDE")

    def test_msa_accepts_a3m_insertions(self):
        stats = validate_a3m(">q\nACDE\n>s\nAgCD-\n", "ACDE")
        self.assertEqual(stats["unique_aligned_rows"], 2)

    def test_control_msa_changes_only_query(self):
        from Bio import SeqIO
        original = ">q\nACGE\n>homolog\nAiCGE\n>second\nACG-\n"
        changed = replace_msa_query(original, "ACGE", "ACAE")
        before = list(SeqIO.parse(io.StringIO(original), "fasta"))
        after = list(SeqIO.parse(io.StringIO(changed), "fasta"))
        self.assertEqual(str(after[0].seq), "ACAE")
        self.assertEqual([str(r.seq) for r in before[1:]], [str(r.seq) for r in after[1:]])

    def test_control_msa_rejects_length_change(self):
        with self.assertRaisesRegex(ValueError, "equal sequence lengths"):
            replace_msa_query(">q\nACGE\n>h\nACG-\n", "ACGE", "ACGEA")

    def test_calibration_requires_every_context_and_seed(self):
        data = {"subjects": dict.fromkeys(CONTEXTS), "seeds": list(SEEDS), "records": [
            {"subject": name, "context": context, "seed": seed, "id": f"{name}-{context}-{seed}"}
            for name, contexts in CONTEXTS.items() for context in contexts for seed in SEEDS]}
        validate_record_coverage(data)
        data["records"] = [r for r in data["records"] if r["context"] != "native_s3p"]
        with self.assertRaisesRegex(ValueError, "Missing"):
            validate_record_coverage(data)

    def test_calibration_duplicate_records_are_not_replicates(self):
        data = {"subjects": dict.fromkeys(CONTEXTS), "seeds": list(SEEDS), "records": [
            {"subject": name, "context": context, "seed": seed, "id": f"{name}-{context}-{seed}"}
            for name, contexts in CONTEXTS.items() for context in contexts for seed in SEEDS]}
        data["records"].append(data["records"][0])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_record_coverage(data)

    def evidence(self):
        source = self.root / "reference.txt"
        source.write_text("Experimental source fixture")
        data = {"target_agi": self.entry.agi,
                "sources": [{"id": "reference", "url": "https://example.org/reference", "artifact": source.name,
                             "sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "scope": "homolog control"}],
                "claims": [{"text": "A homolog control, not target validation.", "source_ids": ["reference"]}]}
        write(self.root / "evidence.json", data)
        return self.root / "evidence.json", data

    def test_cited_evidence_is_not_faked_llm_output(self):
        path, _ = self.evidence()
        result = CitedEvidenceAgent(path).synthesize(self.entry)
        self.assertEqual(result["synthesis"], "curated_evidence_without_llm")
        self.assertIn("[reference]", result["facts"][0])

    def test_evidence_rejects_unknown_citation(self):
        path, _ = self.evidence()
        transport = lambda *_: {"claims": [{"text": "unsupported", "source_ids": ["invented"]}]}
        with self.assertRaisesRegex(ValueError, "cite"):
            CitedEvidenceAgent(path, transport).synthesize(self.entry)

    def test_evidence_rejects_changed_source(self):
        path, _ = self.evidence()
        (self.root / "reference.txt").write_text("tampered")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            CitedEvidenceAgent(path).synthesize(self.entry)

    def quality(self):
        model = self.root / "model.cif"
        model.write_text("coordinate fixture")
        data = {"target_agi": self.entry.agi,
                "sequence_sha256": hashlib.sha256(self.request.target.sequence.encode()).hexdigest(),
                "checks": {k: True for k in ("msa_validated", "chemistry_reviewed", "wt_repeatability",
                                             "experimental_accuracy", "functional_controls", "no_forced_templates")},
                "artifacts": [{"path": model.name, "sha256": hashlib.sha256(model.read_bytes()).hexdigest()}],
                "pocket": {"matched_context": True, "supported_herbicide_positions": [2]}}
        write(self.root / "quality.json", data)
        structures = [replace(self.fixture.structures[0], artifact_path=str(model))]
        return StructureQualityAgent(self.root / "quality.json"), data, structures

    def test_quality_requires_all_checks(self):
        agent, data, structures = self.quality()
        self.assertEqual(agent.assess(self.request, structures)["status"], "PASSED")
        del data["checks"]["experimental_accuracy"]
        write(agent.path, data)
        self.assertEqual(agent.assess(self.request, structures)["status"], "NEEDS_REASSESSMENT")

    def test_quality_binds_to_target_and_artifacts(self):
        agent, _, structures = self.quality()
        with self.assertRaisesRegex(ValueError, "target/sequence"):
            agent.assess(replace(self.request, target=replace(self.request.target, sequence="WRONG")), structures)
        (self.root / "model.cif").write_text("tampered")
        with self.assertRaisesRegex(ValueError, "changed"):
            agent.assess(self.request, structures)

    def test_quality_rejects_uncalibrated_models(self):
        agent, _, structures = self.quality()
        structures = [replace(structures[0], artifact_path=str(self.root / "other.cif"))]
        self.assertEqual(agent.assess(self.request, structures)["status"], "NEEDS_REASSESSMENT")

    def test_pocket_disagreement_blocks(self):
        agent, data, _ = self.quality()
        fingerprint = replace(self.fixture.fingerprints[0], herbicide_selective_mutable=[999])
        self.assertEqual(agent.assess_pocket(self.request, fingerprint, [], [])["status"], "NEEDS_REASSESSMENT")
        data["pocket"]["matched_context"] = False
        write(agent.path, data)
        fingerprint = replace(fingerprint, herbicide_selective_mutable=[2])
        self.assertEqual(agent.assess_pocket(self.request, fingerprint, [], [])["status"], "NEEDS_REASSESSMENT")

    def test_live_workflow_requires_gates(self):
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 5), require_live_gates=True)
        with self.assertRaises(WorkflowReassessmentRequired):
            workflow.run(self.request)

    def test_quality_allows_diagnostics_but_stops_mutation(self):
        quality = Mock()
        quality.assess.return_value = {"status": "NEEDS_REASSESSMENT", "reasons": ["WT variable"]}
        quality.assess_pocket.return_value = {"status": "PASSED", "reasons": []}
        evidence = Mock()
        evidence.synthesize.return_value = {"facts": []}
        dock = Mock(wraps=self.backend)
        workflow = WorkflowOrchestrator(self.registry, self.backend, dock, self.backend, self.backend, None,
                                        ArtifactStore(self.root / "artifacts"), evidence_agent=evidence,
                                        quality_agent=quality, require_live_gates=True)
        workflow.mutation_agent = Mock()
        with self.assertRaisesRegex(WorkflowReassessmentRequired, "WT variable"):
            workflow.run(self.request)
        self.assertTrue(dock.dock.called)
        workflow.mutation_agent.propose.assert_not_called()
        quality.assess_pocket.assert_called_once()
        self.assertTrue(list((self.root / "artifacts").rglob("workflow_stages.json")))

    def test_pocket_stops_before_mutation(self):
        quality = Mock()
        quality.assess.return_value = {"status": "PASSED", "reasons": []}
        quality.assess_pocket.return_value = {"status": "NEEDS_REASSESSMENT", "reasons": ["disputed pocket"]}
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 5), quality_agent=quality)
        workflow.mutation_agent = Mock()
        with self.assertRaisesRegex(WorkflowReassessmentRequired, "disputed pocket"):
            workflow.run(self.request)
        workflow.mutation_agent.propose.assert_not_called()

    def test_reviewer_and_judge_use_separate_backends(self):
        class ReviewOnly(MockScientificBackend):
            def judge(self, *args):
                raise AssertionError("Reviewer used as judge")
        judge = MockScientificBackend()
        judge.model_id = "test-cheap-judge"
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 4), ReviewOnly(), judge_backend=judge)
        self.assertTrue(workflow.run(self.request).judge_results)

    def test_live_judge_does_not_fall_back_to_reviewer(self):
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 5), require_live_gates=True)
        self.assertIsNone(workflow.domain_judge)

    def test_live_packet_has_actual_evidence_and_no_mock_provenance(self):
        old = self.fixture.packets[0]
        prov = [Provenance("local://model", "boltz", "computed")]
        candidate = replace(old.candidate, provenance=prov, metadata={"oracle_evidence": {"predicted_pIC50": 3.2}})
        scores = replace(old.scores, provenance=prov)
        packet = CandidateReviewAgent(None).review(candidate, scores, ["Cited fact"])
        self.assertFalse(any(p.evidence_type == "synthetic" for p in packet.provenance))
        self.assertIn("predicted_pIC50", " ".join(packet.predictions))
        self.assertIn("unavailable", " ".join(packet.unresolved_uncertainty))
        self.assertEqual(packet.status, "NEEDS_REVIEW")

    def test_pose_confidence_not_used_as_affinity(self):
        affinity = Mock()
        affinity.relative_score.return_value = None
        packet = MultiOracleScoringAgent(affinity).score(self.fixture.packets[0].candidate, [], [], self.fixture.fingerprints[0])
        self.assertIn("pose confidence is not affinity", " ".join(packet.uncertainty))

    def test_oracles_require_matching_protocol_and_complete_artifacts(self):
        agent = ArtifactCandidateOracles(self.root / "oracles.json", "protocol")
        with self.assertRaises(WorkflowReassessmentRequired):
            agent.evaluate_candidates(self.request, [self.fixture.packets[0].candidate])
        write(agent.path, {"target_agi": self.entry.agi, "protocol_sha256": "wrong"})
        with self.assertRaisesRegex(ValueError, "protocol mismatch"):
            agent.evaluate_candidates(self.request, [])

    def test_missing_retention_oracles_are_not_filled_from_binding_scores(self):
        agent = ArtifactCandidateOracles(self.root / "absent.json", "protocol")
        self.assertEqual(agent.reference(self.request.target, self.request.herbicide, self.request.native_ligands, []), {})
        self.assertEqual(agent.evaluate(self.request.target, self.fixture.packets[0].candidate,
                                        self.request.herbicide, self.request.native_ligands, []), {})

    def test_api_requires_explicit_model_and_key(self):
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaisesRegex(ValueError, "model ID"):
                OpenAIJSONTransport("")
            with self.assertRaisesRegex(ValueError, "OPENAI_API_KEY"):
                OpenAIJSONTransport("requested-model")

    def test_responses_transport_uses_requested_model_and_no_storage(self):
        captured = []
        def opener(request, **kwargs):
            captured.append(json.loads(request.data))
            return io.BytesIO(json.dumps({"status": "completed", "output": [{"type": "message", "content": [
                {"type": "output_text", "text": json.dumps({"assumptions": [], "uncertainty": ["test"], "recommendation": "more_computation"})}]}]}).encode())
        transport = OpenAIJSONTransport("requested-model", "test-not-a-secret", opener=opener)
        result = transport("candidate_review", {"system_prompt": "test", "facts": []})
        self.assertEqual(result["recommendation"], "more_computation")
        self.assertEqual(captured[0]["model"], "requested-model")
        self.assertFalse(captured[0]["store"])
        self.assertTrue(captured[0]["text"]["format"]["strict"])

    def test_api_incomplete_and_refusal_are_errors(self):
        for response in ({"status": "incomplete"}, {"status": "completed", "output": [{"type": "message", "content": [{"type": "refusal"}]}]}):
            transport = OpenAIJSONTransport("test", "test", opener=lambda *a, **k: io.BytesIO(json.dumps(response).encode()))
            with self.assertRaises(RuntimeError):
                transport("candidate_review", {})

    def test_api_errors_do_not_leak_response_body(self):
        def opener(*args, **kwargs):
            raise HTTPError("https://api.openai.com", 403, "private", {}, io.BytesIO(b"private response"))
        with self.assertRaisesRegex(RuntimeError, "HTTP 403") as caught:
            OpenAIJSONTransport("test", "secret", opener=opener)("candidate_review", {})
        self.assertNotIn("private", str(caught.exception))
        self.assertNotIn("secret", str(caught.exception))

    def test_all_operation_schemas_are_strict(self):
        for operation in ("candidate_review", "domain_judge", "evidence_synthesis"):
            schema = response_schema(operation)
            self.assertFalse(schema["additionalProperties"])
            self.assertEqual(set(schema["required"]), set(schema["properties"]))

    def test_reviewer_receives_retention_failures_before_judging(self):
        reasoner = MockScientificBackend()
        review = Mock(wraps=reasoner.review)
        reasoner.review = review
        workflow = WorkflowOrchestrator(self.registry, *([self.backend] * 4), reasoner)
        workflow.run(self.request)
        predictions = review.call_args_list[0].args[1]
        self.assertIn("INSUFFICIENT_EVIDENCE", " ".join(predictions))

    def test_diagnostic_only_results_cannot_establish_matched_context(self):
        agent, _, _ = self.quality()
        fingerprint = replace(self.fixture.fingerprints[0], herbicide_selective_mutable=[2])
        result = agent.assess_pocket(self.request, fingerprint, [], [], diagnostics={"matched_context": False})
        self.assertEqual(result["status"], "NEEDS_REASSESSMENT")


@unittest.skipUnless(importlib.util.find_spec("rdkit"), "Install calibration extra for pose chemistry tests")
class PoseDiagnosticTests(unittest.TestCase):
    def setUp(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        receptor = self.root / "receptor.pdb"
        receptor.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n")
        self.records = [{"id": "test", "ligand": "test-ligand", "smiles": "CCO", "receptor": receptor.name,
                         "receptor_sha256": hashlib.sha256(receptor.read_bytes()).hexdigest(), "source": str(receptor),
                         "source_sha256": hashlib.sha256(receptor.read_bytes()).hexdigest()}]
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        mol = Chem.RemoveHs(mol)
        self.poses = self.root / "diffdock/test"
        self.poses.mkdir(parents=True)
        for rank in range(1, 5):
            with Chem.SDWriter(str(self.poses / f"rank{rank}_confidence-0.25.sdf")) as writer:
                writer.write(mol)
        self.agent = PoseEnsembleAgent(self.root, self.root / "cache")

    def test_raw_negative_scores_and_numbered_contacts_preserved(self):
        results = self.agent.collect(self.records)
        self.assertEqual(len(results), 4)
        self.assertTrue(all(r["raw_confidence"] == -.25 for r in results))
        self.assertTrue(all(r["contacts"] == [77] for r in results))

    def test_missing_pose_rejected(self):
        (self.poses / "rank2_confidence-0.25.sdf").unlink()
        with self.assertRaisesRegex(ValueError, "four distinct"):
            self.agent.collect(self.records)

    def test_wrong_ligand_rejected(self):
        self.records[0]["smiles"] = "CCN"
        with self.assertRaisesRegex(ValueError, "chemistry mismatch"):
            self.agent.collect(self.records)

    def test_changed_receptor_rejected(self):
        (self.root / "receptor.pdb").write_text("changed")
        with self.assertRaisesRegex(ValueError, "receptor changed"):
            self.agent.collect(self.records)


if __name__ == "__main__":
    unittest.main()
