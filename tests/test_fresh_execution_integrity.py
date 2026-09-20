"""Synthetic I/O and execution-status fixtures; no scientific claims are tested."""
import copy
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from herbicide_desensitization_agent.app.backends.calibration import sha, write
from herbicide_desensitization_agent.app.orchestrator.execution import ExecutionJournal, validate_journal
from herbicide_desensitization_agent.examples.epsps_fresh_campaign import (
    validate_campaign_inputs, analyze, verify_originating_model_hashes, validate_ligand_composition,
)
from herbicide_desensitization_agent.examples.epsps_evidence_stages import literature, merge_journals


class FreshInputIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "msa").mkdir()
        (self.root / "inputs").mkdir()
        self.msa = self.root / "msa/wt.a3m"
        self.msa.write_text(">q\nACDE\n>h\nACD-\n")
        request = {"sequences": [{"protein": {"id": "A", "sequence": "ACDE", "msa": "/remote/campaign/msa/wt.a3m"}}]}
        write(self.root / "inputs/one.yaml", request)
        self.manifest = {"registered_at": "2026-09-01T00:00:00Z", "sequence": "ACDE", "seeds": [1], "protocol": {"version": "test"},
                         "subjects": {"wt": {"sequence": "ACDE", "msa": {"path": "msa/wt.a3m", "sha256": sha(self.msa)}}},
                         "records": [{"id": "one", "subject": "wt", "context": "apo", "mode": "structure", "seed": 1,
                                      "request": "inputs/one.yaml", "request_sha256": sha(self.root / "inputs/one.yaml"), "ligand": None, "ligands": []}]}

    def test_registered_msa_bytes_must_match_even_when_yaml_unchanged(self):
        with patch("herbicide_desensitization_agent.examples.epsps_fresh_campaign.validate_a3m"):
            self.assertEqual(validate_campaign_inputs(self.root, self.manifest)["record_count"], 1)
            self.msa.write_text(">q\nACDE\n>h\nCCDE\n")
            with self.assertRaisesRegex(ValueError, "MSA bytes changed"):
                validate_campaign_inputs(self.root, self.manifest)

    def test_request_identity_is_checked_in_addition_to_hash(self):
        request = json.loads((self.root / "inputs/one.yaml").read_text())
        request["sequences"][0]["protein"]["sequence"] = "AAAA"
        write(self.root / "inputs/one.yaml", request)
        self.manifest["records"][0]["request_sha256"] = sha(self.root / "inputs/one.yaml")
        with patch("herbicide_desensitization_agent.examples.epsps_fresh_campaign.validate_a3m"):
            with self.assertRaisesRegex(ValueError, "sequence or MSA"):
                validate_campaign_inputs(self.root, self.manifest)

    def test_postflight_has_new_measured_event_with_all_source_hashes(self):
        self.manifest["records"][0].update(context="herbicide", mode="affinity", ligand="glyphosate")
        write(self.root / "calibration_inputs.json", self.manifest)
        directory = self.root / "models/one"
        directory.mkdir(parents=True)
        for filename in ("one_model_0.cif", "confidence_one_model_0.json", "affinity_one.json"):
            (directory / filename).write_text("{}")
        effective = copy.deepcopy(self.manifest)
        effective["records"][0]["structure"] = "models/one/one_model_0.cif"
        write(self.root / "effective_model_inputs.json", effective)
        origin = ExecutionJournal(self.root / "model_execution_journal.json")
        origin.execute("affinity", lambda: {"synthetic_fixture": True}, output_artifacts=[directory / "one_model_0.cif"])
        origin.finalize()
        def fake_parser(root, manifest):
            result = [{"id": "one", "synthetic_test_only": True}]
            write(root / "predictions.json", result)
            return result
        with patch("herbicide_desensitization_agent.examples.epsps_fresh_campaign.validate_a3m"), patch(
                "herbicide_desensitization_agent.examples.epsps_fresh_campaign._analyze_predictions", side_effect=fake_parser):
            analyze(self.root)
        journal = json.loads((self.root / "postflight_execution_journal.json").read_text())
        validate_journal(journal)
        event = next(row for row in journal["stages"] if row["stage"] == "model_output_validation")
        names = {Path(row["path"]).name for row in event["input_artifacts"]}
        self.assertTrue({"wt.a3m", "one.yaml", "one_model_0.cif", "confidence_one_model_0.json", "affinity_one.json"} <= names)
        self.assertEqual(event["execution_status"], "COMPLETED")
        self.assertTrue(event["outputs_sha256"])
        self.assertIsNotNone(event["runtime_seconds"])
        self.assertFalse(any(row["stage"] in ("structure", "affinity") for row in journal["stages"]))

    def test_cif_bytes_must_match_unique_originating_event(self):
        directory = self.root / "models/one"
        directory.mkdir(parents=True)
        cif = directory / "one.cif"
        cif.write_text("synthetic bytes\n")
        self.manifest["records"][0]["structure"] = "models/one/one.cif"
        origin_path = self.root / "model_execution_journal.json"
        origin = ExecutionJournal(origin_path)
        origin.execute("structure", lambda: {"synthetic_fixture": True}, output_artifacts=[cif])
        origin.finalize()
        before = origin_path.read_bytes()
        self.assertEqual(verify_originating_model_hashes(self.root, self.manifest)["verified_cif_count"], 1)
        self.assertEqual(origin_path.read_bytes(), before)
        cif.write_text("different downloaded bytes\n")
        with self.assertRaisesRegex(ValueError, "differs from originating"):
            verify_originating_model_hashes(self.root, self.manifest)
        data = json.loads(origin_path.read_text())
        data["stages"][0]["output_artifacts"] = []
        write(origin_path, data)
        with self.assertRaisesRegex(ValueError, "exactly the registered cohort"):
            verify_originating_model_hashes(self.root, self.manifest)

    def test_ambiguous_origin_artifact_is_not_accepted(self):
        directory = self.root / "models/one"
        directory.mkdir(parents=True)
        cif = directory / "one.cif"
        cif.write_text("synthetic bytes\n")
        self.manifest["records"][0]["structure"] = "models/one/one.cif"
        second = copy.deepcopy(self.manifest["records"][0])
        second.update(id="two", structure="models/two/two.cif")
        self.manifest["records"].append(second)
        origin = ExecutionJournal(self.root / "model_execution_journal.json")
        origin.execute("structure", lambda: {}, output_artifacts=[cif, cif])
        origin.finalize()
        with self.assertRaisesRegex(ValueError, "unambiguous"):
            verify_originating_model_hashes(self.root, self.manifest)


class LigandInventoryTests(unittest.TestCase):
    @staticmethod
    def structure(chains):
        class Model:
            def __init__(self):
                self.chains = {name: types.SimpleNamespace(id=name, get_atoms=lambda symbols=symbols: iter([types.SimpleNamespace(element=e) for e in symbols])) for name, symbols in chains.items()}
            def __iter__(self):
                return iter(self.chains.values())
            def __getitem__(self, name):
                return self.chains[name]
        return [Model()]

    def chemistry(self):
        atoms = [types.SimpleNamespace(GetSymbol=lambda: "C", GetAtomicNum=lambda: 6), types.SimpleNamespace(GetSymbol=lambda: "O", GetAtomicNum=lambda: 8)]
        molecule = types.SimpleNamespace(GetAtoms=lambda: atoms, GetNumHeavyAtoms=lambda: 2)
        return patch.dict(sys.modules, {"rdkit": types.SimpleNamespace(Chem=types.SimpleNamespace(MolFromSmiles=lambda text: molecule))})

    def test_missing_or_wrong_ligand_atoms_rejected(self):
        record = {"id": "fixture", "ligands": [{"id": "B", "name": "test", "model_smiles": "CO"}]}
        with self.chemistry():
            result = validate_ligand_composition(self.structure({"A": ["C"], "B": ["C", "O", "H"]}), record)
            self.assertEqual(result["B"]["observed_heavy_atoms"], 2)
            for symbols in (["C"], ["C", "N"], ["C", "O", "O"]):
                with self.assertRaisesRegex(ValueError, "heavy-atom inventory"):
                    validate_ligand_composition(self.structure({"A": ["C"], "B": symbols}), record)

    def test_apo_requires_absence_of_ligand_chains(self):
        with self.chemistry():
            self.assertEqual(validate_ligand_composition(self.structure({"A": ["C"]}), {"id": "apo", "ligands": []}), {})
            with self.assertRaisesRegex(ValueError, "including apo absence"):
                validate_ligand_composition(self.structure({"A": ["C"], "B": ["O"]}), {"id": "apo", "ligands": []})


class EvidenceExecutionTests(unittest.TestCase):
    def test_literature_failure_does_not_suppress_independent_chemistry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sources"
            source.mkdir()
            baseline = root / "baseline.json"
            write(baseline, {"sources": {"test": {"url": "https://example.invalid/test-fixture"}}})
            for name in ("PMC8983318.xml", "Funke09.pdf", "7PXY.cif"):
                (source / name).write_text("synthetic I/O test only")
            write(source / "chemical_states_ph74.json", {"ligands": {"test": {"selected_campaign_state": {"smiles": "C"}}}})
            chem = types.SimpleNamespace(MolFromSmiles=lambda s: s, FindMolChiralCenters=lambda mol, **kw: [],
                                         MolToSmiles=lambda s: s, GetFormalCharge=lambda mol: 0)
            with patch("herbicide_desensitization_agent.examples.epsps_evidence_stages.resolve_api_key", side_effect=RuntimeError("credential unavailable")), patch.dict(sys.modules, {"rdkit": types.SimpleNamespace(Chem=chem)}):
                with self.assertRaisesRegex(RuntimeError, "both attempted"):
                    literature(root / "campaign", baseline, source)
            data = json.loads((root / "campaign/evidence/execution_journal.json").read_text())
            validate_journal(data)
            states = {row["stage"]: row["execution_status"] for row in data["stages"]}
            self.assertEqual(states, {"literature": "FAILED", "chemical_state": "COMPLETED"})
            self.assertEqual(data["status"], "COMPLETED_WITH_FAILURES")

    def test_merge_preserves_actual_timing_aliases_and_dependency_identities(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = ExecutionJournal(root / "one.json", run_id="first")
            first.execute("structure", lambda: {"actual_fixture_call": True})
            first.execute("stability_scoring", lambda: {"actual_fixture_call": True}, dependencies=["structure"])
            original = copy.deepcopy(first.latest("stability_scoring"))
            first.finalize()
            second = ExecutionJournal(root / "two.json", run_id="second")
            second.skip("mutation_design", "SKIPPED_NO_CANDIDATES", "No accepted candidates")
            second.finalize()
            result = merge_journals([root / "one.json", root / "two.json"], root / "merged.json")
            self.assertEqual(result["status"], "COMPLETED_WITH_SKIPS")
            stability = next(row for row in result["stages"] if row["stage"] == "stability")
            self.assertEqual(stability["source_stage"], "stability_scoring")
            self.assertEqual(stability["source_event_id"], original["id"])
            for key in ("started_at", "finished_at", "runtime_seconds", "inputs_sha256", "outputs_sha256"):
                self.assertEqual(stability[key], original[key])
            ids = {row["id"] for row in result["stages"]}
            self.assertTrue(set(stability["dependencies"]) <= ids)
            self.assertEqual(stability["source_dependencies"], ["structure"])


if __name__ == "__main__":
    unittest.main()
