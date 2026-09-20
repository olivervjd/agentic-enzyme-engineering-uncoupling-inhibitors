"""Preregistered thresholds derived before candidate inspection, never fitted to candidates."""
from datetime import datetime
from .schema import digest, finite, scientific_provenance

WEAKENING_DIRECTIONS = {"Kd": "increase", "Ki": "increase", "IC50": "increase",
                        "pIC50": "decrease", "pKd": "decrease", "pKi": "decrease", "delta_G": "increase"}


def calibration_errors(calibration):
    if not isinstance(calibration, dict):
        return ["calibration missing"]
    errors = []
    payload = {k: v for k, v in calibration.items() if k != "digest"}
    if calibration.get("digest") != digest(payload):
        errors.append("calibration digest mismatch (modified after freezing)")
    if calibration.get("status") != "VALIDATED":
        errors.append("WT/control calibration not validated")
    if not calibration.get("thresholds") or not calibration.get("wt_variability"):
        errors.append("thresholds or WT variability missing")
    try:
        recomputed = calibrate_thresholds(calibration["specification"], calibration["wt_replicates"],
                                         calibration["input_controls"], registered_at=calibration["registered_at"],
                                         candidates_started_at=calibration["candidates_started_at"])
        if recomputed != calibration:
            errors.append("calibration does not reproduce from registered WT/control inputs")
    except (KeyError, ValueError, TypeError):
        errors.append("calibration source data missing or invalid")
    return errors


def calibrate_thresholds(spec, wt_replicates, controls, *, registered_at, candidates_started_at):
    """spec={metrics:{name:{direction:min|max,limit,margin,unit}}, binding:{...}}.

    The empirical WT envelope plus preregistered margin is clipped to the
    preregistered acceptable limit. Controls have expected=pass|fail, kind,
    metrics and provenance. Any required control missing or misclassified
    prevents validation. A new calibration requires a new prospective run.
    """
    before = datetime.fromisoformat(registered_at.replace("Z", "+00:00"))
    after = datetime.fromisoformat(candidates_started_at.replace("Z", "+00:00"))
    if before.tzinfo is None or after.tzinfo is None or before >= after:
        raise ValueError("Calibration must be timezone-aware and preregistered before candidate evaluation")
    if not spec.get("metrics"):
        raise ValueError("Preregistered metric specifications required")
    for name, rule in spec["metrics"].items():
        if rule.get("direction") not in ("min", "max") or not finite(rule.get("limit")) or not finite(rule.get("margin")) or rule["margin"] < 0 or not rule.get("unit"):
            raise ValueError(f"Invalid preregistered threshold: {name}")
        if name in ("tm_score_forward", "tm_score_reverse", "ca_lddt", "coverage") and rule["direction"] != "min":
            raise ValueError(f"Structural similarity requires a lower acceptance bound: {name}")
        if ("rmsd" in name or name in ("folding_ddg", "pocket_volume_change", "contact_change", "assembly_change")) and rule["direction"] != "max":
            raise ValueError(f"Structural perturbation requires an upper acceptance bound: {name}")
    errors, thresholds, variability = [], {}, {}
    binding = spec.get("binding", {})
    required_evaluations = spec.get("required_evaluations")
    if not isinstance(required_evaluations, dict) or not required_evaluations:
        errors.append("authoritative preregistered evaluation registry missing")
        required_evaluations = {}
    for key in ("minimum_replicates", "minimum_methods"):
        if not isinstance(binding.get(key), int) or binding[key] < 2:
            errors.append(f"{key} must be preregistered and at least two")
    for key in ("pose_consistency_min", "contact_reproducibility_min"):
        if not finite(binding.get(key)) or not 0 < binding[key] <= 1:
            errors.append(f"{key} must be a preregistered probability")
    for evaluation_id, role in required_evaluations.items():
        rule = binding.get("evaluations", {}).get(evaluation_id, {})
        if not rule or not rule.get("unit") or rule.get("metric") not in WEAKENING_DIRECTIONS:
            errors.append(f"ligand-specific affinity rule missing: {evaluation_id}")
            continue
        values = []
        for row in wt_replicates:
            measurement = row.get("binding_measurements", {}).get(evaluation_id, {})
            if measurement.get("metric") != rule["metric"] or measurement.get("unit") != rule["unit"] or not finite(measurement.get("value")):
                errors.append(f"WT ligand replicate measurement missing: {evaluation_id}")
                break
            values.append(measurement["value"])
        if len(values) == len(wt_replicates) and values:
            observed_range = max(values)-min(values)
            variability[evaluation_id] = {"min": min(values), "max": max(values), "range": observed_range,
                                         "unit": rule["unit"], "n": len(values)}
            if not finite(rule.get("wt_variability_bound")) or rule["wt_variability_bound"] < observed_range:
                errors.append(f"declared WT variability understates observed range: {evaluation_id}")
        if role == "herbicide" and (rule.get("weakening_direction") != WEAKENING_DIRECTIONS[rule["metric"]]
                                    or not finite(rule.get("weakening_threshold")) or rule["weakening_threshold"] <= 0):
            errors.append(f"incorrect physical weakening direction/threshold: {evaluation_id}")
    if "herbicide" not in required_evaluations.values() or "native_substrate" not in required_evaluations.values():
        errors.append("registered herbicide and native substrate evaluations required")
    seen = set()
    for row in wt_replicates:
        key = (row.get("method"), row.get("seed"))
        if None in key or key in seen or not scientific_provenance(row.get("provenance")):
            errors.append("WT replicates need unique method/seed and scientific provenance")
        seen.add(key)
    minimum_wt_replicates = spec.get("minimum_wt_replicates", 3)
    if not isinstance(minimum_wt_replicates, int) or minimum_wt_replicates < 3 or len(wt_replicates) < minimum_wt_replicates:
        errors.append("insufficient independent WT replicates")
    for name, rule in spec["metrics"].items():
        values = [row.get("metrics", {}).get(name) for row in wt_replicates]
        if not values or not all(finite(x) for x in values):
            errors.append(f"missing WT metric: {name}")
            continue
        if name in ("pocket_volume_change", "contact_change", "assembly_change"):
            values = [abs(value) for value in values]
        variability[name] = {"min": min(values), "max": max(values), "range": max(values)-min(values),
                             "n": len(values), "unit": rule["unit"]}
        bound = min(rule["limit"], max(values)+rule["margin"]) if rule["direction"] == "max" else max(rule["limit"], min(values)-rule["margin"])
        thresholds[name] = {"direction": rule["direction"], "limit": bound, "unit": rule["unit"]}
        if any((v > bound if rule["direction"] == "max" else v < bound) for v in values):
            errors.append(f"WT baseline fails preregistered acceptance for {name}")
    rows = []
    required_kinds = {"resistant", "catalytic_loss", "positive_binding", "negative_binding"}
    required_kinds.update(spec.get("required_control_kinds", []))
    if required_kinds - {r.get("kind") for r in controls}:
        errors.append("missing required control classes: " + ", ".join(sorted(required_kinds-{r.get("kind") for r in controls})))
    for row in controls:
        metrics = {name: abs(value) if name in ("pocket_volume_change", "contact_change", "assembly_change") and finite(value) else value
                   for name, value in row.get("metrics", {}).items()}
        complete = (bool(thresholds) and all(finite(metrics.get(n)) for n in thresholds)
                    and scientific_provenance(row.get("provenance")))
        passed = complete and all((metrics[n] <= r["limit"] if r["direction"] == "max" else metrics[n] >= r["limit"]) for n, r in thresholds.items())
        predicted = ("pass" if passed else "fail") if complete else "unknown"
        # Class labels alone never validate a panel. Check ligand-resolved
        # measurements recover the documented functional control phenotype.
        ligand_checks = []
        kind = row.get("kind")
        for evaluation_id, role in required_evaluations.items():
            rule = binding.get("evaluations", {}).get(evaluation_id, {})
            if kind in ("resistant", "catalytic_loss", "neutral"):
                measurement = row.get("ligand_changes", {}).get(evaluation_id, {})
                valid = (measurement.get("metric") == rule.get("metric") and measurement.get("unit") == rule.get("unit")
                         and finite(measurement.get("value")))
                if not valid:
                    ligand_checks.append(False)
                    continue
                value = measurement["value"]
                if role == "herbicide" and kind == "resistant":
                    bound = rule.get("weakening_threshold")
                    ligand_checks.append(finite(bound) and (value > bound if rule.get("weakening_direction") == "increase" else value < -bound))
                elif role == "native_substrate":
                    interval = rule.get("equivalence_interval")
                    valid = isinstance(interval, list) and len(interval) == 2 and all(finite(x) for x in interval)
                    equivalent = valid and interval[0] <= value <= interval[1]
                    ligand_checks.append(valid and (not equivalent if kind == "catalytic_loss" else equivalent))
            elif kind in ("positive_binding", "negative_binding"):
                measurement = row.get("ligand_affinities", {}).get(evaluation_id, {})
                interval = rule.get("binding_support_interval")
                valid = (measurement.get("metric") == rule.get("metric") and measurement.get("unit") == rule.get("unit")
                         and finite(measurement.get("value")) and isinstance(interval, list) and len(interval) == 2
                         and all(finite(x) for x in interval))
                if valid:
                    supports = interval[0] <= measurement["value"] <= interval[1]
                    weakened = (measurement["value"] > interval[1] if WEAKENING_DIRECTIONS.get(rule.get("metric")) == "increase"
                                else measurement["value"] < interval[0])
                    ligand_checks.append(supports if kind == "positive_binding" else weakened)
                else:
                    ligand_checks.append(False)
        ligand_correct = bool(ligand_checks) and all(ligand_checks)
        correct = complete and row.get("expected") == predicted and ligand_correct
        rows.append({**row, "predicted": predicted, "ligand_controls_recovered": ligand_correct, "correct": correct})
        if not correct:
            errors.append(f"control not recovered: {row.get('id', 'unnamed')}")
    result = {"schema_version": "1.0", "status": "VALIDATED" if not errors else "INSUFFICIENT_EVIDENCE",
              "registered_at": registered_at, "candidates_started_at": candidates_started_at,
              "specification": spec, "thresholds": thresholds, "wt_variability": variability,
              "wt_replicates": wt_replicates, "input_controls": controls,
              "controls": rows, "limitations": errors,
              "binding": spec.get("binding", {}), "contact": spec.get("contact", {})}
    result["digest"] = digest(result)
    return result


def validate_live_registry(calibration, *, target_id, herbicide, native_ligands, epsps=False):
    """Authoritative request-to-frozen-protocol binding before mutation generation."""
    errors = calibration_errors(calibration)
    spec = (calibration or {}).get("specification", {})
    if spec.get("target_id") != target_id or spec.get("herbicide_id") != herbicide:
        errors.append("calibration target/herbicide differs from the live request")
    registered_native = set(spec.get("native_ligand_ids", []))
    if set(native_ligands) != registered_native:
        errors.append("all registered native substrates must appear in calibration scope")
    required_baseline = ("experimental_structure_comparison", "published_binding_site", "catalytic_residues",
                         "documented_msa", "biological_assembly", "repeated_structure_prediction", "repeated_docking")
    baseline = spec.get("baseline_evidence", {})
    for name in required_baseline:
        entry = baseline.get(name, {})
        if entry.get("validated") is not True or not scientific_provenance(entry.get("provenance")):
            errors.append(f"WT baseline evidence incomplete: {name}")
    neutral = baseline.get("neutral_mutation_controls", {})
    if not (neutral.get("validated") is True and scientific_provenance(neutral.get("provenance"))) and not (neutral.get("unavailable") is True and neutral.get("search_provenance")):
        errors.append("neutral mutation controls require evidence or a documented unavailable search")
    if epsps:
        represented = set(spec.get("evaluated_complexes", []))
        required = {"glyphosate+S3P", "PEP+S3P", "S3P+PEP", "apo", "substrate_bound_control"}
        if not required <= represented:
            errors.append("EPSPS requires glyphosate+S3P, PEP+S3P, S3P+PEP, apo and substrate-bound controls")
    return sorted(set(errors))
