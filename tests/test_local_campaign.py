import unittest
import json
import tempfile
from pathlib import Path

from herbicide_desensitization_agent.examples.epsps_local_campaign import CONDITIONS, docking_contact_support, mutate, report, validate_affinity, write_json


class LocalCampaignTests(unittest.TestCase):
    def test_wt_and_exact_single_substitution(self):
        sequence = "A" * 520
        self.assertEqual(mutate(sequence, "WT"), sequence)
        result = mutate(sequence, "A288V")
        self.assertEqual(result[287], "V")
        self.assertEqual(sum(a != b for a, b in zip(sequence, result)), 1)

    def test_mutations_outside_domain_or_wrong_identity_fail(self):
        for mutation in ("A76V", "M288L", "A521V", "A288A"):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                mutate("A" * 520, mutation)

    def test_native_substrate_affinities_keep_same_complex_context(self):
        self.assertEqual(CONDITIONS["native"], ("pep", "B"))
        self.assertEqual(CONDITIONS["native_s3p"], ("pep", "C"))
        self.assertEqual(CONDITIONS["herbicide"], ("glyphosate", "B"))

    def test_affinity_units_and_required_values(self):
        self.assertEqual(validate_affinity({"affinity_pred_value": 2., "affinity_probability_binary": .4}), 4.)
        with self.assertRaises(ValueError):
            validate_affinity({"affinity_pred_value": 2.})

    def test_invalid_affinity_outputs_are_rejected(self):
        for invalid in (float("nan"), float("inf"), True, "2"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                validate_affinity({"affinity_pred_value": invalid, "affinity_probability_binary": .4})
        for invalid in (-.1, 1.1):
            with self.subTest(probability=invalid), self.assertRaises(ValueError):
                validate_affinity({"affinity_pred_value": 2., "affinity_probability_binary": invalid})

    def test_failed_run_cannot_publish_success_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_json(root / "execution_manifest.json", {"status": "FAILED"})
            with self.assertRaisesRegex(ValueError, "incomplete or failed"):
                report(root)
            self.assertFalse((root / "CAMPAIGN_REPORT.md").exists())

    def test_docking_contact_support_uses_wt_not_mutant_poses(self):
        poses = [{"mutation": "WT", "ligand": "glyphosate", "contacts_uniprot": [288]},
                 {"mutation": "WT", "ligand": "pep", "contacts_uniprot": []},
                 {"mutation": "M288L", "ligand": "pep", "contacts_uniprot": [288]}]
        row = docking_contact_support(poses, ["WT", "M288L", "M288I"])[0]
        self.assertEqual(row, {"sequence_position": 288, "glyphosate_poses": 1, "glyphosate_contacts": 1,
                               "pep_poses": 1, "pep_contacts": 0})

    def test_report_preserves_wt_ranges_and_blocked_overall_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_json(root / "execution_manifest.json", {"status": "COMPUTATION_COMPLETE_ANALYSIS_PENDING"})
            write_json(root / "inputs/input_manifest.json", {"mutations": ["WT"]})
            write_json(root / "affinity_predictions.json", {"records": [
                {"mutation": "WT", "ligand": ligand, "replicate": replicate, "predicted_pIC50": value}
                for ligand in ("glyphosate", "pep", "s3p") for replicate, value in ((1, 1.), (2, 3.))]})
            write_json(root / "docking_predictions.json", {"records": []})
            assessment = root / "assessment/AT2G45300-glyphosate"
            write_json(assessment / "retention_evidence.json", {"WT": {"structural_metrics": {
                "query_tm_score": 1., "target_tm_score": 1., "alignment_lddt": 1., "alignment_coverage": 1.,
                "global_ca_rmsd_angstrom": 0., "active_site_rmsd_angstrom": 0.}}})
            write_json(assessment / "function_retention_report.json", [{"mutation": "WT", "decision": "REFERENCE"}])
            report(root, git_commit="test-fixture")
            text = (root / "CAMPAIGN_REPORT.md").read_text()
            self.assertIn("2.000 [1.000, 3.000]", text)
            self.assertIn("not confidence intervals", text)
            state = json.loads((root / "execution_manifest.json").read_text())
            self.assertEqual(state["overall_status"], "INCOMPLETE_BLOCKED_STAGES")
