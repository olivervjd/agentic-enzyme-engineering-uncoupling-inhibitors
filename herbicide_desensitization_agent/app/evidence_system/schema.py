"""JSON contracts and deterministic validators for evidence-aware screening.

Unknown values are null, never zero. Scientific values must carry tool or
experimental provenance; language-model explanations are not measurements.
"""
import hashlib
import json
import math

EVIDENCE_STATES = (
    "Experimentally confirmed binding",
    "Experimentally confirmed non-binding or substantially weakened binding",
    "Strong computational support for binding",
    "Strong computational support for weakened binding",
    "Conflicting computational evidence", "Insufficient evidence", "Not evaluated",
)
DECISIONS = (
    "MEETS_COMPUTATIONAL_SCREEN", "FAILS_COMPUTATIONAL_SCREEN",
    "INSUFFICIENT_EVIDENCE", "METHOD_DISAGREEMENT_REQUIRES_REVIEW",
    "EXPERIMENTALLY_VALIDATED", "EXPERIMENTALLY_REJECTED",
)
RESIDUE_CLASSES = (
    "HERBICIDE_SELECTIVE_CONTACT", "SHARED_HERBICIDE_NATIVE_CONTACT",
    "NATIVE_CRITICAL_CONTACT", "CATALYTIC_OR_PROTECTED", "SECOND_SHELL_CANDIDATE",
    "UNSTABLE_OR_METHOD_DEPENDENT_CONTACT", "INSUFFICIENT_EVIDENCE",
)
CONTEXT_FIELDS = (
    "protein_sequence", "isoform", "modeled_residue_range", "residue_mapping",
    "oligomeric_state", "cellular_context", "substrates", "cofactors", "ligands",
    "pH", "ionic_assumptions", "metal_ions", "catalytic_waters", "partners",
    "structure_model", "structure_model_version", "affinity_model",
    "affinity_model_version", "msa_source", "msa_depth", "seed", "sampling",
)
CHEMICAL_FIELDS = ("id", "canonical_identifier", "stereochemistry", "protonation",
                   "tautomer", "formal_charge")
STRUCTURAL_METRICS = ("tm_score_forward", "tm_score_reverse", "ca_lddt", "coverage",
                      "global_ca_rmsd", "active_site_backbone_rmsd", "pocket_sidechain_rmsd",
                      "ligand_rmsd", "pocket_volume_change", "contact_change", "folding_ddg")


def finite(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def scientific_provenance(value):
    """Traceable imported computational/experimental values only; no LLM numerics."""
    return (isinstance(value, dict) and value.get("source_type") in ("computational", "experimental")
            and bool(value.get("source")) and bool(value.get("artifact")))


def context_errors(context):
    if not isinstance(context, dict):
        return ["context"]
    errors = [key for key in CONTEXT_FIELDS if key not in context or context[key] is None]
    unknown = {"", "unknown", "unspecified", "unavailable", "pending", "missing"}
    for key in ("protein_sequence", "isoform", "oligomeric_state", "cellular_context", "ionic_assumptions",
                "structure_model", "structure_model_version", "affinity_model", "affinity_model_version", "msa_source"):
        value = context.get(key)
        if not isinstance(value, str) or value.strip().lower() in unknown:
            errors.append(f"resolved {key}")
    for key in ("substrates", "cofactors", "ligands", "metal_ions", "catalytic_waters", "partners"):
        if not isinstance(context.get(key), list):
            errors.append(f"explicit list {key}")
    if not isinstance(context.get("residue_mapping"), dict) or not context.get("residue_mapping"):
        errors.append("explicit residue_mapping")
    modeled_range = context.get("modeled_residue_range")
    if not isinstance(modeled_range, list) or len(modeled_range) != 2 or not all(isinstance(x, int) and not isinstance(x, bool) for x in modeled_range) or not 1 <= modeled_range[0] <= modeled_range[1]:
        errors.append("ordered positive modeled_residue_range")
    if not isinstance(context.get("sampling"), dict) or not context.get("sampling"):
        errors.append("explicit sampling settings")
    for i, ligand in enumerate(context.get("ligands") if isinstance(context.get("ligands"), list) else []):
        if not isinstance(ligand, dict):
            errors.append(f"ligands[{i}]")
            continue
        errors.extend(f"ligands[{i}].{key}" for key in CHEMICAL_FIELDS
                      if key not in ligand or ligand[key] is None)
        errors.extend(f"ligands[{i}].resolved {key}" for key in CHEMICAL_FIELDS if key != "formal_charge"
                      and (not isinstance(ligand.get(key), str) or ligand[key].strip().lower() in unknown))
        if not isinstance(ligand.get("formal_charge"), int) or isinstance(ligand.get("formal_charge"), bool):
            errors.append(f"ligands[{i}].integer formal_charge")
    if not context.get("ligands"):
        # Apo is represented explicitly, not by an unknown/missing ligand field.
        if context.get("complex_type") != "apo":
            errors.append("ligands or explicit complex_type=apo")
    if "pH" in context and (not finite(context["pH"]) or not 0 <= context["pH"] <= 14):
        errors.append("finite pH")
    for key in ("msa_depth", "seed"):
        if key in context and (not isinstance(context[key], int) or isinstance(context[key], bool) or context[key] < (1 if key == "msa_depth" else 0)):
            errors.append(f"integer {key}")
    return sorted(set(errors))


def compare_contexts(left, right, *, mutant=False, cross_method=False):
    """No inference of chemical equivalence. Lists for sets are canonicalized.

    Cross-method comparisons permit different model identities and seeds, but
    require them to be recorded. Paired WT/mutant comparisons require identity
    of models, seeds and sampling, and explicit canonical substitutions.
    """
    missing = [f"{side}.{x}" for side, ctx in (("left", left), ("right", right))
               for x in context_errors(ctx)]
    left, right = left if isinstance(left, dict) else {}, right if isinstance(right, dict) else {}
    ignore = {"protein_sequence"} if mutant else set()
    if cross_method:
        ignore.update(("structure_model", "structure_model_version", "affinity_model",
                       "affinity_model_version", "seed"))
    set_fields = {"substrates", "cofactors", "ligands", "metal_ions", "catalytic_waters", "partners"}
    def normalized(key, value):
        if key in set_fields and isinstance(value, list):
            return sorted(json.dumps(x, sort_keys=True) for x in value)
        return value
    mismatches = [key for key in CONTEXT_FIELDS if key not in ignore
                  and normalized(key, left.get(key)) != normalized(key, right.get(key))]
    if mutant:
        a, b = left.get("protein_sequence", ""), right.get("protein_sequence", "")
        changes = right.get("substitutions")
        actual = [f"{aa}{i}{bb}" for i, (aa, bb) in enumerate(zip(a, b), 1) if aa != bb]
        if len(a) != len(b) or not changes or sorted(changes) != sorted(actual):
            mismatches.append("declared canonical substitutions")
    return {"equivalent": not missing and not mismatches, "missing": missing,
            "mismatches": mismatches, "comparison": "cross_method" if cross_method else "matched_pair"}


def affinity_errors(affinity):
    if not isinstance(affinity, dict):
        return ["affinity missing"]
    errors = []
    allowed = {"Kd": {"M", "mM", "uM", "nM", "pM"}, "Ki": {"M", "mM", "uM", "nM", "pM"},
               "IC50": {"M", "mM", "uM", "nM", "pM"}, "pIC50": {"dimensionless"},
               "pKd": {"dimensionless"}, "pKi": {"dimensionless"}, "delta_G": {"kcal/mol", "kJ/mol"}}
    if affinity.get("metric") not in allowed or affinity.get("unit") not in allowed.get(affinity.get("metric"), set()):
        errors.append("unsupported physical metric/unit; no confidence-to-affinity conversion")
    interval = affinity.get("interval")
    if not isinstance(interval, (list, tuple)) or len(interval) != 2 or not all(finite(x) for x in interval):
        errors.append("finite uncertainty interval missing")
    elif not finite(affinity.get("estimate")) or not interval[0] <= affinity["estimate"] <= interval[1]:
        errors.append("estimate outside ordered interval")
    elif affinity.get("metric") in ("Kd", "Ki", "IC50") and interval[0] <= 0:
        errors.append("concentration interval must be positive")
    if not affinity.get("interval_method") or not finite(affinity.get("interval_level")) or not 0 < affinity.get("interval_level", 0) < 1:
        errors.append("uncertainty method/confidence level missing")
    if not scientific_provenance(affinity.get("provenance")):
        errors.append("traceable scientific affinity provenance missing")
    return errors


STAGE_REQUIRED = {
    "registry": ("target", "required_evaluations", "required_context"),
    "literature": ("structures", "binding_site", "catalytic_residues", "controls", "references"),
    "chemical_state": ("context",),
    "structure": ("replicates", "metrics", "provenance"),
    "docking": ("replicates", "pose_clusters", "provenance"),
    "affinity": ("affinity",), "contacts": ("contacts", "definitions", "provenance"),
    "calibration": ("thresholds", "wt_variability", "controls", "digest"),
    "mutation_design": ("mutations", "rationales"), "stability": ("folding_ddg", "provenance"),
    "function_retention": ("decision", "gates"), "evidence_review": ("claims", "limitations"),
    "dashboard": ("tables", "manifest"), "learning_loop": ("assays", "review_status"),
}


def validate_stage(stage, payload):
    if stage not in STAGE_REQUIRED:
        raise ValueError(f"Unknown stage: {stage}")
    if not isinstance(payload, dict):
        return {"stage": stage, "valid": False, "errors": ["object payload required"]}
    errors = [f"missing {x}" for x in STAGE_REQUIRED[stage] if x not in payload or payload[x] is None]
    if stage == "chemical_state":
        errors.extend(context_errors(payload.get("context")))
    if stage == "affinity":
        errors.extend(affinity_errors(payload.get("affinity")))
    if "provenance" in STAGE_REQUIRED[stage] and not scientific_provenance(payload.get("provenance")):
        errors.append("scientific provenance required")
    if stage in ("structure", "docking"):
        rows = payload.get("replicates", [])
        keys = {(r.get("method"), r.get("seed")) for r in rows if isinstance(r, dict) and r.get("method") and r.get("seed") is not None}
        if len(keys) < 2 or len(keys) != len(rows):
            errors.append("at least two unique method/seed replicates required")
        if stage == "structure" and not payload.get("metrics"):
            errors.append("nonempty structural measurements required")
        if stage == "docking" and not payload.get("pose_clusters"):
            errors.append("pose clustering evidence required")
    if stage == "calibration":
        from .calibration import calibration_errors
        errors.extend(calibration_errors(payload))
    if stage == "registry" and any(not payload.get(key) for key in STAGE_REQUIRED[stage]):
        errors.append("nonempty target, required evaluations and biological context registry required")
    if stage == "function_retention" and (payload.get("decision") not in DECISIONS or not payload.get("gates")):
        errors.append("recognized decision and explicit evidence gates required")
    return {"stage": stage, "valid": not errors, "errors": sorted(set(errors))}
