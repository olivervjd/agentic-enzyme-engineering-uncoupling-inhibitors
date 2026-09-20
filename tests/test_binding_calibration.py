"""Synthetic values in this file test software only, never scientific performance."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from herbicide_desensitization_agent.app.backends.binding_calibration import (
    fit_calibration, calibrated_prediction, control_diagnostics, load_dataset, build_binding_report,
)
from herbicide_desensitization_agent.app.backends.credentials import resolve_api_key


def fixture():
    scope = dict(target_id="TEST_ONLY", assay_id="fixture", conditions={"pH": 7},
                 prediction_protocol="fixture", ligand_scope="test series")
    rows = []
    for split, n in (("train", 8), ("calibration", 9), ("test", 10)):
        for i in range(n):
            x = 2 + i * 4 / (n - 1)
            y = x + .6 + (.1 if split == "calibration" else 0)
            ident = f"{split}-{i}"
            rows.append(dict(id=ident, subject=ident, ligand="fixture", group=ident, split=split,
                             scope=scope.copy(), endpoint="IC50", relation="=", unit="M", value=10 ** -y,
                             predicted_pIC50=x, evidence_type="experimental"))
    return dict(scope=scope, observations=rows)


class BindingCalibrationTests(unittest.TestCase):
    def test_affine_fit_and_heldout_interval(self):
        report = fit_calibration(fixture())
        self.assertEqual(report["status"], "VALIDATED_WITHIN_SCOPE")
        self.assertAlmostEqual(report["model"]["intercept"], .6)
        self.assertAlmostEqual(report["model"]["interval_half_width"], .1)
        self.assertEqual(report["metrics"]["test_interval_coverage"], 1)
        self.assertFalse(report["functional_retention_established"])
        self.assertAlmostEqual(calibrated_prediction(report, 3, report["scope"])["estimate"], 3.6)

    def test_does_not_fit_on_test(self):
        data = fixture()
        before = fit_calibration(data)
        for r in data["observations"]:
            if r["split"] == "test":
                r["value"] *= 100
        after = fit_calibration(data)
        self.assertEqual(before["model"], after["model"])
        self.assertEqual(after["status"], "FAILED_VALIDATION")

    def test_incompatible_endpoints_and_censoring_excluded(self):
        for endpoint, relation in (("Kd", "="), ("Ki", "="), ("Km", "="), ("IC50", ">")):
            data = fixture()
            for r in data["observations"]:
                r.update(endpoint=endpoint, relation=relation)
            report = fit_calibration(data)
            self.assertEqual(report["eligible"], 0)
            self.assertIsNone(report["model"])

    def test_scope_mismatch_excluded(self):
        data = fixture()
        data["observations"][0]["scope"] = {"target_id": "wrong"}
        self.assertEqual(fit_calibration(data)["status"], "INSUFFICIENT_DATA")

    def test_dependent_replicates_rejected(self):
        for key in ("id", "group", "subject"):
            data = fixture()
            data["observations"][8][key] = data["observations"][0][key]
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                fit_calibration(data)

    def test_missing_scope_and_small_panel_abstain(self):
        self.assertEqual(fit_calibration({})["status"], "INSUFFICIENT_DATA")
        data = fixture()
        data["observations"] = data["observations"][:2]
        self.assertIsNone(fit_calibration(data)["model"])

    def test_rejects_nan_zero_and_infinity(self):
        for value in (0, -1, float("nan"), float("inf"), True):
            data = fixture()
            data["observations"][0]["value"] = value
            with self.assertRaises(ValueError):
                fit_calibration(data)

    def test_does_not_extrapolate_or_transfer_scope(self):
        report = fit_calibration(fixture())
        with self.assertRaisesRegex(ValueError, "Extrapolation"):
            calibrated_prediction(report, 9, report["scope"])
        with self.assertRaisesRegex(ValueError, "scope"):
            calibrated_prediction(report, 3, {})

    def test_zero_dynamic_range(self):
        data = fixture()
        for row in data["observations"]:
            row["predicted_pIC50"] = 3
        self.assertIsNone(fit_calibration(data)["model"])

    def test_controls_are_diagnostics_not_numeric_calibration(self):
        rows = [dict(subject=subject, ligand=ligand, seed=seed, predicted_pIC50=value)
                for subject, ligand, value in (("ecoli_wt", "glyphosate", 4), ("ecoli_G96A", "glyphosate", 3),
                                               ("ecoli_wt", "pep", 4), ("ecoli_G96A", "pep", 5))
                for seed in (101, 103, 107)]
        results = control_diagnostics(rows)
        self.assertEqual(results[0]["status"], "DIRECTION_CONSISTENT")
        self.assertEqual(results[1]["status"], "DIRECTION_DISCORDANT")
        self.assertFalse(results[0]["numerically_comparable"])
        self.assertEqual(control_diagnostics([])[0]["status"], "MISSING_MATCHED_CONTROLS")

    def test_source_traversal_rejected(self):
        from herbicide_desensitization_agent.app.backends.binding_calibration import _artifact
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "inside"):
                _artifact(Path(tmp), {"path": "../secret", "sha256": "x"})

    def test_units_are_normalized(self):
        for unit, scale in (("mM", 1000), ("uM", 1e6), ("nM", 1e9)):
            data = fixture()
            for row in data["observations"]:
                row.update(unit=unit, value=row["value"] * scale)
            self.assertAlmostEqual(fit_calibration(data)["model"]["intercept"], .6)

    def test_one_subject_cannot_cross_splits_under_different_ligands(self):
        data = fixture()
        data["observations"][8].update(subject=data["observations"][0]["subject"], ligand="different")
        with self.assertRaisesRegex(ValueError, "leak"):
            fit_calibration(data)

    def test_control_duplicate_seed_rejected(self):
        row = dict(subject="ecoli_wt", ligand="glyphosate", seed=101, predicted_pIC50=4)
        with self.assertRaisesRegex(ValueError, "Duplicate control"):
            control_diagnostics([row, row])

    def test_loader_derives_prediction_and_checks_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.txt"
            source.write_text("Synthetic test only")
            predicted = root / "predictions.json"
            predicted.write_text(json.dumps({"records": [dict(id="p", subject="test", ligand="test", predicted_pIC50=4)]}))
            dataset = fixture()
            row = dataset["observations"][0]
            row.update(subject="test", ligand="test", predicted_pIC50=99,
                       source={"url": "https://example.org", "locator": "Table 1", "path": "source.txt",
                               "sha256": hashlib.sha256(source.read_bytes()).hexdigest()},
                       prediction={"path": "predictions.json", "sha256": hashlib.sha256(predicted.read_bytes()).hexdigest(), "record_ids": ["p"]})
            dataset["observations"] = [row]
            path = root / "dataset.json"
            path.write_text(json.dumps(dataset))
            self.assertEqual(load_dataset(path)["observations"][0]["predicted_pIC50"], 4)
            source.write_text("Changed bytes")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                load_dataset(path)

    def test_first_calibration_run_can_verify_in_memory_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "predictions.json").write_text('{"records": []}')
            (root / "calibration_inputs.json").write_text('{"target_agi":"AT2G45300"}')
            manifest = {"artifacts": [{"path": n, "sha256": hashlib.sha256((root / n).read_bytes()).hexdigest()}
                                      for n in ("predictions.json", "calibration_inputs.json")]}
            report = build_binding_report(root, input_report=manifest)
            self.assertEqual(report["status"], "INSUFFICIENT_DATA")
            self.assertEqual(report["controls"][0]["status"], "CONFOUNDED_CONTROL_PROTOCOL")
            (root / "predictions.json").write_text('{"records": [1]}')
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                build_binding_report(root, input_report=manifest)

    def test_negative_slope_fails_validation(self):
        data = fixture()
        for row in data["observations"]:
            row["value"] = 10 ** -(8-row["predicted_pIC50"])
        report = fit_calibration(data)
        self.assertEqual(report["status"], "FAILED_VALIDATION")
        self.assertFalse(report["checks"]["positive_slope"])

    @patch.dict("os.environ", {}, clear=True)
    def test_codex_api_key_reuse_requires_opt_in_and_never_oauth(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "auth.json"
            path.write_text(json.dumps({"auth_mode": "apikey", "OPENAI_API_KEY": "test-key"}))
            self.assertIsNone(resolve_api_key(codex_home=tmp)[0])
            self.assertEqual(resolve_api_key(use_codex=True, codex_home=tmp), ("test-key", "codex_api_key"))
            path.write_text(json.dumps({"auth_mode": "chatgpt", "tokens": {"access_token": "never-use"}}))
            self.assertIsNone(resolve_api_key(use_codex=True, codex_home=tmp)[0])

    @patch.dict("os.environ", {"OPENAI_API_KEY": "env-test"})
    def test_environment_key_preferred(self):
        self.assertEqual(resolve_api_key(use_codex=True), ("env-test", "environment"))


if __name__ == "__main__":
    unittest.main()
