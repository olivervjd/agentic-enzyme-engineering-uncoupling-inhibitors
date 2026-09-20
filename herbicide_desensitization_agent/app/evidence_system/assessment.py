"""Deterministic ligand-specific evidence ladder and conservative screen decisions."""
from .schema import (EVIDENCE_STATES, STRUCTURAL_METRICS, affinity_errors, compare_contexts,
                     finite, scientific_provenance)
from .calibration import calibration_errors, WEAKENING_DIRECTIONS


def assess_binding(record, calibration=None):
    """record: ligand, variant, methods[], optional experimental evidence.

    Each method records independent_family, context, conclusion, replicates
    [{seed,artifact}], affinity, pose_consistency, contact_reproducibility and
    control_comparison. Replicate count counts unique method-family/seed pairs.
    """
    result = {"variant": record.get("variant"), "ligand": record.get("ligand"),
              "evidence_state": EVIDENCE_STATES[5], "method_agreement": "not_assessed",
              "independent_replicates": 0, "limitations": [], "methods": record.get("methods", [])}
    experimental = record.get("experimental")
    if experimental and scientific_provenance(experimental.get("provenance")) and experimental.get("provenance", {}).get("source_type") == "experimental" and experimental.get("context_matched") is True:
        measurement = experimental.get("measurement")
        direct_measurement = (not affinity_errors(measurement) and measurement.get("metric") in ("Kd", "pKd")
                              and measurement["provenance"]["source_type"] == "experimental"
                              and experimental.get("direct_binding_assay") is True)
        comparison = experimental.get("comparison", {})
        if experimental.get("outcome") == "weakened_or_nonbinding":
            # A negative call additionally requires a quantitative matched control.
            delta = affinity_difference(comparison.get("wt_affinity"), measurement)
            threshold = comparison.get("threshold")
            direct_measurement = direct_measurement and delta["valid"] and finite(threshold) and threshold > 0
            direct_measurement = direct_measurement and comparison.get("wt_affinity", {}).get("provenance", {}).get("source_type") == "experimental"
            if direct_measurement:
                direct_measurement = delta["interval"][0] > threshold if measurement["metric"] == "Kd" else delta["interval"][1] < -threshold
        if experimental.get("outcome") in ("binding", "weakened_or_nonbinding") and experimental.get("assay") and direct_measurement:
            result["evidence_state"] = EVIDENCE_STATES[0 if experimental["outcome"] == "binding" else 1]
            result["experimental"] = experimental
            return result
    methods = record.get("methods", [])
    if not methods:
        result["evidence_state"] = EVIDENCE_STATES[6] if not record.get("evaluated") else EVIDENCE_STATES[5]
        result["limitations"] = ["No usable ligand-specific evidence"]
        return result
    issues = calibration_errors(calibration)
    binding_rules = (calibration or {}).get("binding", {})
    for key in ("minimum_replicates", "minimum_methods", "pose_consistency_min", "contact_reproducibility_min"):
        if key not in binding_rules:
            issues.append(f"preregistered {key} missing")
    families, replicate_keys, conclusions = set(), set(), set()
    equivalent = True
    for i, method in enumerate(methods):
        family = method.get("independent_family")
        if not family:
            issues.append(f"method {i}: independence family missing")
        else:
            families.add(family)
        comparison = compare_contexts(methods[0].get("context"), method.get("context"), cross_method=True)
        if not comparison["equivalent"]:
            equivalent = False
            issues.append(f"method {i}: non-equivalent or incomplete chemical/biological context")
        issues.extend(f"method {i}: {x}" for x in affinity_errors(method.get("affinity")))
        keys = set()
        artifacts = set()
        for rep in method.get("replicates", []):
            if rep.get("seed") is None or not rep.get("artifact"):
                issues.append(f"method {i}: replicate provenance missing")
            else:
                keys.add((family, rep["seed"]))
                artifacts.add(rep["artifact"])
        if len(keys) != len(method.get("replicates", [])) or len(artifacts) != len(keys):
            issues.append(f"method {i}: duplicate or invalid replicates")
        replicate_keys.update(keys)
        minimum_replicates = binding_rules.get("minimum_replicates", 2)
        if not isinstance(minimum_replicates, int) or len(keys) < max(2, minimum_replicates):
            issues.append(f"method {i}: insufficient independent replicates")
        for name in ("pose_consistency", "contact_reproducibility"):
            value, limit = method.get(name), binding_rules.get(name + "_min")
            if not finite(value) or not 0 <= value <= 1 or not finite(limit) or value < limit:
                issues.append(f"method {i}: {name} not established")
        control = method.get("control_comparison", {})
        if control.get("passed") is not True or not scientific_provenance(control.get("provenance")) or not control.get("positive_control_id") or not control.get("negative_control_id"):
            issues.append(f"method {i}: positive/negative control comparison missing or failed")
        rule = binding_rules.get("evaluations", {}).get(record.get("evaluation_id"), {})
        affinity = method.get("affinity", {})
        if (rule.get("metric"), rule.get("unit")) != (affinity.get("metric"), affinity.get("unit")):
            issues.append(f"method {i}: ligand-specific calibrated endpoint missing or mismatched")
        elif method.get("conclusion") == "supports_weakening":
            delta = affinity_difference(method.get("wt_affinity"), affinity)
            context = compare_contexts(method.get("wt_context"), method.get("context"), mutant=True)
            bound, noise = rule.get("weakening_threshold"), rule.get("wt_variability_bound")
            if not delta["valid"] or not context["equivalent"] or not finite(bound) or bound <= 0 or not finite(noise):
                issues.append(f"method {i}: calibrated WT-relative weakening evidence missing")
            else:
                supported = (delta["interval"][0] > max(bound, noise) if WEAKENING_DIRECTIONS.get(affinity.get("metric")) == "increase"
                             else delta["interval"][1] < -max(bound, noise))
                if not supported:
                    issues.append(f"method {i}: weakening unsupported throughout uncertainty interval")
        elif method.get("conclusion") == "supports_binding":
            interval, support = affinity.get("interval"), rule.get("binding_support_interval")
            if (not isinstance(support, list) or len(support) != 2 or not all(finite(x) for x in support)
                    or affinity_errors(affinity) or interval[0] < support[0] or interval[1] > support[1]):
                issues.append(f"method {i}: affinity interval outside or missing calibrated binding-support interval")
        if method.get("conclusion") in ("supports_binding", "supports_weakening"):
            conclusions.add(method["conclusion"])
        else:
            issues.append(f"method {i}: uncertain conclusion")
    result["independent_replicates"] = len(replicate_keys)
    minimum_methods = binding_rules.get("minimum_methods", 2)
    if not isinstance(minimum_methods, int) or len(families) < max(2, minimum_methods):
        issues.append("insufficient independent methods")
    result["method_agreement"] = "non_equivalent" if not equivalent else "disagreement" if len(conclusions) > 1 else "agreement" if len(families) >= 2 and len(conclusions) == 1 else "insufficient_methods"
    if equivalent and len(conclusions) > 1:
        result["evidence_state"] = EVIDENCE_STATES[4]
    elif not issues:
        result["evidence_state"] = EVIDENCE_STATES[2 if conclusions == {"supports_binding"} else 3]
    result["limitations"] = sorted(set(issues))
    return result


def affinity_difference(wt, mutant):
    """Conservative interval subtraction; not a fitted paired confidence interval.

    No transformation between physical quantities or units is performed.
    """
    errors = affinity_errors(wt) + affinity_errors(mutant)
    if wt and mutant and (wt.get("metric"), wt.get("unit")) != (mutant.get("metric"), mutant.get("unit")):
        errors.append("WT/mutant physical quantities or units differ")
    if errors:
        return {"valid": False, "limitations": sorted(set(errors))}
    return {"valid": True, "metric": wt["metric"], "unit": wt["unit"],
            "estimate": mutant["estimate"]-wt["estimate"],
            "interval": [mutant["interval"][0]-wt["interval"][1], mutant["interval"][1]-wt["interval"][0]],
            "interval_method": "conservative endpoint subtraction; joint coverage not asserted",
            "limitations": ["Marginal interval subtraction does not establish paired statistical confidence"]}


def evaluate_candidate(candidate, calibration=None):
    gates = []
    def gate(name, status, reason):
        gates.append({"name": name, "status": status, "reason": reason})
    errors = calibration_errors(calibration)
    gate("wt_calibration", "MISSING" if errors else "PASS", "; ".join(errors) if errors else "Frozen WT/control calibration validated")
    required = candidate.get("required_evaluations")
    registered = (calibration or {}).get("specification", {}).get("required_evaluations", {})
    if not registered or not required or set(required) != set(registered):
        gate("registered_evaluations", "MISSING", "Required evaluations must match the frozen target/native-function registry")
    records = candidate.get("bindings", [])
    names = [x.get("evaluation_id") for x in records]
    if not required or len(names) != len(set(names)):
        gate("evidence_completeness", "MISSING", "Explicit required ligand/complex registry and unique evaluation IDs required")
    else:
        missing = set(required)-set(names)
        gate("evidence_completeness", "MISSING" if missing else "PASS", "Missing: " + ", ".join(sorted(missing)) if missing else "All registered complexes present")
    roles = {r.get("role") for r in records}
    if "herbicide" not in roles or "native_substrate" not in roles:
        gate("required_roles", "MISSING", "Both herbicide and native-substrate evaluations required")
    binding_rules = (calibration or {}).get("binding", {}).get("evaluations", {})
    for row in records:
        key = row.get("evaluation_id", "unknown")
        if row.get("role") != registered.get(key):
            gate(f"{key}:registered_role", "MISSING", "Evaluation role absent or differs from frozen registry")
        context = compare_contexts(row.get("wt_context"), row.get("mutant_context"), mutant=True)
        gate(f"{key}:matched_context", "PASS" if context["equivalent"] else "MISSING", "Matched WT/mutant context" if context["equivalent"] else str(context))
        assessment = assess_binding(row.get("assessment", {}), calibration)
        source_assessment = row.get("assessment", {})
        coupled = (source_assessment.get("variant") == candidate.get("variant")
                   and source_assessment.get("ligand") == row.get("ligand")
                   and source_assessment.get("evaluation_id") == key)
        matching_methods = [m for m in source_assessment.get("methods", [])
                            if m.get("affinity") == row.get("mutant_affinity") and m.get("wt_affinity") == row.get("wt_affinity")
                            and compare_contexts(row.get("mutant_context"), m.get("context"))["equivalent"]
                            and compare_contexts(row.get("wt_context"), m.get("wt_context"))["equivalent"]]
        coupled = coupled and bool(matching_methods)
        for method in source_assessment.get("methods", []):
            coupled = coupled and compare_contexts(row.get("mutant_context"), method.get("context"), cross_method=True)["equivalent"]
        if not coupled:
            gate(f"{key}:evidence_identity", "MISSING", "Assessment must belong to this variant, ligand, complex, context and imported affinity")
        if assessment["method_agreement"] == "disagreement":
            gate(f"{key}:method_agreement", "DISAGREEMENT", "Independent methods conflict under equivalent context")
        elif assessment["evidence_state"] not in EVIDENCE_STATES[:4]:
            gate(f"{key}:binding_evidence", "MISSING", "; ".join(assessment["limitations"]) or assessment["evidence_state"])
        else:
            gate(f"{key}:binding_evidence", "PASS", assessment["evidence_state"])
        delta = affinity_difference(row.get("wt_affinity"), row.get("mutant_affinity"))
        rule = binding_rules.get(key)
        if not delta["valid"] or not rule or (rule.get("metric"), rule.get("unit")) != (delta.get("metric"), delta.get("unit")):
            gate(f"{key}:affinity", "MISSING", "Matched affinity intervals and preregistered ligand-specific rule required")
            continue
        lo, hi = delta["interval"]
        if row.get("role") == "herbicide":
            direction, bound = rule.get("weakening_direction"), rule.get("weakening_threshold")
            noise = rule.get("wt_variability_bound")
            if direction != WEAKENING_DIRECTIONS.get(delta["metric"]) or not finite(bound) or bound <= 0 or not finite(noise) or noise < 0:
                gate(f"{key}:weakening", "MISSING", "Calibrated weakening direction, positive threshold and WT noise required")
            else:
                passes = lo > max(bound, noise) if direction == "increase" else hi < -max(bound, noise)
                gate(f"{key}:weakening", "PASS" if passes else "FAIL", f"Entire change interval {delta['interval']} must exceed threshold {bound} and WT variability {noise} {delta['unit']}")
        else:
            interval = rule.get("equivalence_interval")
            if not isinstance(interval, list) or len(interval) != 2 or not all(finite(x) for x in interval) or interval[0] > interval[1]:
                gate(f"{key}:native_equivalence", "MISSING", "Preregistered native equivalence interval required")
            else:
                gate(f"{key}:native_equivalence", "PASS" if lo >= interval[0] and hi <= interval[1] else "FAIL", f"Entire change interval {delta['interval']} inside {interval} {delta['unit']}")
    metrics = candidate.get("structural_metrics", {})
    thresholds = (calibration or {}).get("thresholds", {})
    assembly_required = (calibration or {}).get("specification", {}).get("assembly_required")
    if not isinstance(assembly_required, bool):
        gate("assembly_requirement", "MISSING", "Frozen registry must state whether assembly evaluation is required")
    for name in STRUCTURAL_METRICS + (("assembly_change",) if assembly_required else ()):
        value, rule = metrics.get(name), thresholds.get(name)
        if not isinstance(value, dict) or not finite(value.get("value")) or not scientific_provenance(value.get("provenance")) or not rule or value.get("unit") != rule.get("unit"):
            gate(name, "MISSING", "Scientific metric, matching unit and calibrated threshold required")
        else:
            metric = value["value"]
            if ("rmsd" in name and metric < 0) or (name in ("tm_score_forward", "tm_score_reverse", "coverage") and not 0 <= metric <= 1):
                gate(name, "MISSING", "Metric outside physical domain")
                continue
            if name in ("pocket_volume_change", "contact_change", "assembly_change"):
                metric = abs(metric)
            passes = metric <= rule["limit"] if rule["direction"] == "max" else metric >= rule["limit"]
            gate(name, "PASS" if passes else "FAIL", f"{value['value']} {rule['direction']} bound {rule['limit']} {rule['unit']}")
    for name in ("biological_context_complete", "contacts_reproducible", "protected_interactions_retained"):
        evidence = candidate.get(name)
        if not isinstance(evidence, dict) or evidence.get("passed") not in (True, False) or not scientific_provenance(evidence.get("provenance")):
            gate(name, "MISSING", "Traceable tool/experimental gate evidence required")
        else:
            gate(name, "PASS" if evidence["passed"] else "FAIL", evidence.get("reason", "Evidence recorded"))
    states = {g["status"] for g in gates}
    decision = ("METHOD_DISAGREEMENT_REQUIRES_REVIEW" if "DISAGREEMENT" in states else
                "FAILS_COMPUTATIONAL_SCREEN" if "FAIL" in states else
                "INSUFFICIENT_EVIDENCE" if "MISSING" in states else "MEETS_COMPUTATIONAL_SCREEN")
    assay = candidate.get("experimental_validation")
    if assay and assay.get("scope") == "herbicide_interference_and_native_function" and assay.get("reviewed") is True and scientific_provenance(assay.get("provenance")) and assay["provenance"]["source_type"] == "experimental" and assay.get("measurement_artifacts"):
        if assay.get("outcome") in ("validated", "rejected"):
            decision = "EXPERIMENTALLY_VALIDATED" if assay["outcome"] == "validated" else "EXPERIMENTALLY_REJECTED"
    return {"variant": candidate.get("variant"), "decision": decision, "gates": gates,
            "limitations": [g["reason"] for g in gates if g["status"] != "PASS"],
            "claim": "Computational screen only; herbicide resistance and retained enzyme function require biochemical validation"}
