"""Endpoint-matched affine calibration with independent split-conformal residuals.

This estimates IC50 only. It cannot establish Kd, native turnover, or resistance.
"""
from __future__ import annotations

import hashlib
import json
import math
import statistics as stats
from pathlib import Path


LEGEND = (
    "Boltz predicted pIC50 = 6 - affinity_pred_value; larger values mean stronger predicted inhibition. "
    "Seed ranges show computational variability, not experimental confidence intervals. IC50, Ki, Km and "
    "Kd are not interchangeable. The numeric calibrator uses measured IC50 only, in a single target, "
    "assay-condition and prediction-protocol scope. A linear fit uses training groups, absolute residuals "
    "from separate calibration groups define a 90% split-conformal prediction interval, and untouched test "
    "groups assess error and coverage. Coverage relies on exchangeability within that scope, not arbitrary "
    "new mutations. Passing this diagnostic does not establish substrate retention or enzyme function."
)
POLICY = {"minimum_train": 8, "minimum_calibration": 9, "minimum_test": 10,
          "alpha": .1, "maximum_test_rmse": .5, "minimum_test_coverage": .9,
          "maximum_interval_half_width": 1.0}
SCOPE_FIELDS = ("target_id", "assay_id", "conditions", "prediction_protocol", "ligand_scope")


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Expected a finite numeric measurement")
    return float(value)


def _artifact(root, spec):
    path = (root / spec["path"]).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Calibration source must be inside the dataset directory")
    if hashlib.sha256(path.read_bytes()).hexdigest() != spec["sha256"]:
        raise ValueError("Calibration source hash mismatch")
    return path


def load_dataset(path):
    """Verify evidence bytes and derive each prediction from referenced model outputs."""
    path = Path(path)
    dataset = json.loads(path.read_text())
    for row in dataset["observations"]:
        source = row["source"]
        if not source["url"].startswith("https://") or not source.get("locator"):
            raise ValueError("Experimental measurements need a source URL and table/page locator")
        _artifact(path.parent, source)
        predicted = json.loads(_artifact(path.parent, row["prediction"]).read_text())
        ids = row["prediction"]["record_ids"]
        rows = [r for r in predicted["records"] if r["id"] in ids]
        if not ids or len(set(ids)) != len(ids) or len(rows) != len(ids):
            raise ValueError("Prediction records are absent or duplicated")
        if any(r.get("subject") != row["subject"] or r.get("ligand") != row["ligand"] for r in rows):
            raise ValueError("Prediction subject/ligand mismatch")
        row["predicted_pIC50"] = stats.mean(_number(r["predicted_pIC50"]) for r in rows)
    dataset["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dataset


def fit_calibration(dataset):
    """Predeclared train/calibration/test groups, never randomly split seed replicates."""
    report = {"status": "INSUFFICIENT_DATA", "endpoint": "pIC50", "policy": dict(POLICY),
              "scope": dataset.get("scope"), "input_sha256": dataset.get("input_sha256"),
              "eligible": 0, "excluded": [], "reasons": [], "model": None, "legend": LEGEND,
              "functional_retention_established": False}
    scope = dataset.get("scope", {})
    if any(not scope.get(k) for k in SCOPE_FIELDS):
        report["reasons"].append("A target, assay, conditions, prediction protocol and ligand scope are required")
        return report
    groups, identities, ids = set(), set(), set()
    subject_splits = {}
    splits = {name: [] for name in ("train", "calibration", "test")}
    factors = {"M": 1., "mM": 1e-3, "uM": 1e-6, "nM": 1e-9}
    for row in dataset.get("observations", []):
        identity = (row["subject"], row["ligand"])
        if row["id"] in ids or row["group"] in groups or identity in identities:
            raise ValueError("Duplicate observation or dependent group; aggregate replicates before splitting")
        ids.add(row["id"])
        groups.add(row["group"])
        identities.add(identity)
        if row["subject"] in subject_splits and subject_splits[row["subject"]] != row["split"]:
            raise ValueError("A protein subject cannot leak across data splits")
        subject_splits[row["subject"]] = row["split"]
        reason = None
        if row.get("endpoint") != "IC50":
            reason = "Endpoint mismatch: only measured IC50 is eligible"
        elif row.get("relation") != "=":
            reason = "Censored measurements are not eligible for this point calibrator"
        elif row.get("scope") != scope:
            reason = "Different target, assay conditions or prediction protocol"
        elif row.get("evidence_type") != "experimental":
            reason = "Experimental measurement required"
        if reason:
            report["excluded"].append({"id": row["id"], "reason": reason})
            continue
        if row["split"] not in splits or row["unit"] not in factors:
            raise ValueError("Invalid split or IC50 concentration unit")
        measured = _number(row["value"]) * factors[row["unit"]]
        if measured <= 0:
            raise ValueError("IC50 must be positive")
        splits[row["split"]].append({"id": row["id"], "x": _number(row["predicted_pIC50"]),
                                     "y": -math.log10(measured)})
    report["counts"] = {k: len(v) for k, v in splits.items()}
    report["eligible"] = sum(report["counts"].values())
    for name in splits:
        if len(splits[name]) < POLICY["minimum_" + name]:
            report["reasons"].append(f"Need at least {POLICY['minimum_' + name]} independent {name} groups")
    if report["reasons"]:
        return report
    train, calibration, test = (splits[k] for k in ("train", "calibration", "test"))
    mx, my = (stats.mean(r[k] for r in train) for k in ("x", "y"))
    denominator = sum((r["x"] - mx) ** 2 for r in train)
    if denominator < 1e-8:
        report["reasons"].append("Training predictions have no usable dynamic range")
        return report
    slope = sum((r["x"] - mx) * (r["y"] - my) for r in train) / denominator
    intercept = my - slope * mx
    predict = lambda r: intercept + slope * r["x"]
    residuals = sorted(abs(r["y"] - predict(r)) for r in calibration)
    rank = math.ceil((len(calibration) + 1) * (1 - POLICY["alpha"]))
    width = residuals[rank - 1]
    rmse = math.sqrt(stats.mean((r["y"] - predict(r)) ** 2 for r in test))
    raw_rmse = math.sqrt(stats.mean((r["y"] - r["x"]) ** 2 for r in test))
    coverage = stats.mean(abs(r["y"] - predict(r)) <= width + 1e-12 for r in test)
    training_range = [min(r["x"] for r in train), max(r["x"] for r in train)]
    checks = {"positive_slope": slope > 0, "test_rmse": rmse <= POLICY["maximum_test_rmse"],
              "not_worse_than_raw": rmse <= raw_rmse + 1e-12,
              "test_coverage": coverage >= POLICY["minimum_test_coverage"],
              "interval_precision": width <= POLICY["maximum_interval_half_width"],
              "heldout_in_training_range": all(training_range[0] <= r["x"] <= training_range[1] for r in calibration + test)}
    report.update(status="VALIDATED_WITHIN_SCOPE" if all(checks.values()) else "FAILED_VALIDATION",
                  checks=checks, model={"slope": slope, "intercept": intercept, "interval_half_width": width,
                                       "prediction_range": training_range},
                  metrics={"test_rmse": rmse, "raw_test_rmse": raw_rmse, "test_interval_coverage": coverage},
                  test_points=[{**r, "calibrated": predict(r), "low": predict(r)-width, "high": predict(r)+width} for r in test])
    report["reasons"] = [k + " failed" for k, ok in checks.items() if not ok]
    return report


def calibrated_prediction(report, value, scope):
    value = _number(value)
    if report["status"] != "VALIDATED_WITHIN_SCOPE" or scope != report["scope"]:
        raise ValueError("No validated calibration for this scope")
    model = report["model"]
    if not model["prediction_range"][0] <= value <= model["prediction_range"][1]:
        raise ValueError("Extrapolation beyond calibration range is prohibited")
    estimate = model["intercept"] + model["slope"] * value
    return {"endpoint": "pIC50", "estimate": estimate,
            "interval": [estimate-model["interval_half_width"], estimate+model["interval_half_width"]], "legend": LEGEND}


def control_diagnostics(predictions):
    results = []
    for ligand in ("glyphosate", "pep"):
        values = {}
        for subject in ("ecoli_wt", "ecoli_G96A"):
            rows = [r for r in predictions if r["subject"] == subject and r["ligand"] == ligand]
            if len({r["seed"] for r in rows}) != len(rows):
                raise ValueError("Duplicate control seed")
            values[subject] = {r["seed"]: _number(r["predicted_pIC50"]) for r in rows}
        wt, mutant = values["ecoli_wt"], values["ecoli_G96A"]
        if set(wt) != {101, 103, 107} or set(wt) != set(mutant):
            results.append({"ligand": ligand, "status": "MISSING_MATCHED_CONTROLS"})
            continue
        delta = [mutant[s] - wt[s] for s in sorted(wt)]
        results.append({"ligand": ligand, "status": "DIRECTION_CONSISTENT" if max(delta) < 0 else "DIRECTION_DISCORDANT",
                        "wild_type_mean": stats.mean(wt.values()), "mutant_mean": stats.mean(mutant.values()),
                        "paired_delta_mean": stats.mean(delta), "paired_delta_range": [min(delta), max(delta)],
                        "expected_direction": "weaker in G96A", "numerically_comparable": False,
                        "source": "https://www.rcsb.org/structure/1MI4",
                        "scope": "E. coli homolog qualitative control, not Arabidopsis validation"})
    return results


def build_binding_report(root, dataset_path=None, *, input_report=None):
    root = Path(root)
    prediction_path = root / "predictions.json"
    structural_report = input_report if input_report is not None else json.loads((root / "calibration_report.json").read_text())
    recorded = {a["path"]: a for a in structural_report["artifacts"]}
    for name in ("predictions.json", "calibration_inputs.json"):
        if name not in recorded:
            raise ValueError("Calibration manifest lacks a required binding input")
        _artifact(root, recorded[name])
    target_id = json.loads((root / "calibration_inputs.json").read_text())["target_agi"]
    dataset = load_dataset(dataset_path) if dataset_path else None
    if dataset and (dataset.get("scope", {}).get("target_id") != target_id or
                    dataset.get("application_prediction_sha256") != hashlib.sha256(prediction_path.read_bytes()).hexdigest()):
        raise ValueError("Binding dataset is not bound to this target and prediction ensemble")
    predictions = json.loads(prediction_path.read_text())["records"]
    report = fit_calibration(dataset) if dataset else {
        "status": "INSUFFICIENT_DATA", "endpoint": "pIC50", "eligible": 0, "model": None,
        "reasons": ["No endpoint-matched measured IC50 calibration dataset supplied",
                    "G96A kinetic observations cannot be used as pIC50 labels"],
        "policy": dict(POLICY), "functional_retention_established": False, "legend": LEGEND}
    report["controls"] = control_diagnostics(predictions)
    if structural_report.get("control_msa_homolog_rows_matched") is not True:
        report["controls"] = [{"ligand": c["ligand"], "status": "CONFOUNDED_CONTROL_PROTOCOL"} for c in report["controls"]]
    report["prediction_sha256"] = hashlib.sha256(prediction_path.read_bytes()).hexdigest()
    report["target_id"] = target_id
    report["functional_retention_reason"] = "Requires matched substrate-binding and activity evidence; calibrated IC50 alone cannot pass the function gate"
    return report
