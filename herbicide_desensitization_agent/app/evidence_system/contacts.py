"""Coordinate-derived contacts with explicit geometry, denominators and limitations.

PDB atom positions are measurements from the input file. Chemical atom typing
must come from a validated preparer; untyped PDBs yield proximity/clash/metal
coordination hypotheses only, never invented hydrogen bonds or salt bridges.
"""
from collections import defaultdict
from math import acos, degrees, dist, sqrt
from .schema import finite, scientific_provenance

CONTACT_DEFINITIONS = {
    "general_proximity": {"max_distance_angstrom": 4.5, "definition": "Any non-hydrogen protein/ligand atom pair"},
    "hydrogen_bond": {"max_distance_angstrom": 3.5, "minimum_D_H_A_angle_degrees": 120,
                      "definition": "Typed donor and acceptor with explicit hydrogen geometry"},
    "salt_bridge": {"max_distance_angstrom": 4.0, "definition": "Explicitly typed opposite formal-charge atoms"},
    "hydrophobic_contact": {"max_distance_angstrom": 4.5, "definition": "Both atoms explicitly typed hydrophobic"},
    "metal_coordination": {"max_distance_angstrom": 2.8, "definition": "Metal to N/O/S; coordination hypothesis without electronic validation"},
    "steric_clash": {"max_fraction_vdw_sum": 0.65, "definition": "Distance < 0.65 times sum of tabulated van der Waals radii; excludes typed bonded pairs"},
    "aromatic_interaction": {"definition": "Not assigned without ring plane geometry; missing in PDB-only extraction"},
    "water_mediated": {"definition": "Not assigned without explicit water and both hydrogen-bond geometries"},
}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75}
METALS = {"MG", "MN", "ZN", "FE", "CA", "CO", "NI", "CU"}


def parse_pdb_atoms(text, *, model=1):
    """Keep a single MODEL, highest-occupancy alternate location (ties blank/A).

    PDB fixed-width coordinates only; unsupported/malformed records are errors,
    never silently reassigned. PDB insertion codes are retained in residue IDs.
    """
    atoms, selected_model = {}, 1
    for line in text.splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            selected_model = int(line[10:14].strip())
        if record not in ("ATOM", "HETATM") or selected_model != model:
            continue
        if len(line) < 54:
            raise ValueError("Truncated PDB atom record")
        name, resname = line[12:16].strip(), line[17:20].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = next((c.upper() for c in name if c.isalpha()), "")
        if element in ("H", "D"):
            continue
        xyz = [float(line[a:b]) for a,b in ((30,38), (38,46), (46,54))]
        if not all(finite(x) for x in xyz):
            raise ValueError("Non-finite atom coordinate")
        residue = line[22:26].strip() + line[26:27].strip()
        chain, alt = line[21:22].strip(), line[16:17].strip()
        occupancy = float(line[54:60].strip() or 1) if len(line) >= 60 else 1.0
        atom = {"name": name, "element": element, "xyz": xyz, "chain": chain,
                "structure_residue": residue, "amino_acid": resname, "record": record,
                "occupancy": occupancy, "altloc": alt}
        key = (record, chain, residue, resname, name)
        score = (occupancy, 2 if alt == "" else 1 if alt == "A" else 0)
        if key not in atoms or score > atoms[key][0]:
            atoms[key] = (score, atom)
    return [atoms[k][1] for k in sorted(atoms)]


def _hbond(donor, acceptor):
    if not donor.get("donor") or not acceptor.get("acceptor"):
        return False
    for hydrogen in donor.get("hydrogens", []):
        a = [x-y for x,y in zip(donor["xyz"], hydrogen)]
        b = [x-y for x,y in zip(acceptor["xyz"], hydrogen)]
        denominator = sqrt(sum(x*x for x in a)*sum(x*x for x in b))
        if denominator and degrees(acos(max(-1,min(1,sum(x*y for x,y in zip(a,b))/denominator)))) >= 120:
            return True
    return False


def extract_pose_contacts(protein_atoms, ligand_atoms, *, variant, ligand, pose_id,
                          replicate_id, model_id, method, canonical_mapping,
                          annotations=None, provenance=None):
    """Return one row per contacting residue/ligand atom pair.

    canonical_mapping keys are 'chain:structure_residue', values canonical
    integer residue numbers. Unknown mappings remain null and block nomination.
    Annotations are keyed by canonical residue as str and never inferred.
    """
    if not scientific_provenance(provenance):
        raise ValueError("Coordinate artifact provenance required")
    for atom in [*protein_atoms, *ligand_atoms]:
        if not isinstance(atom.get("xyz"), (list, tuple)) or len(atom["xyz"]) != 3 or not all(finite(x) for x in atom["xyz"]):
            raise ValueError("Finite Cartesian coordinates required")
    rows = []
    for atom in protein_atoms:
        if atom.get("element") in ("H", "D"):
            continue
        canonical = canonical_mapping.get(f"{atom['chain']}:{atom['structure_residue']}")
        annotation = (annotations or {}).get(str(canonical), {})
        for other in ligand_atoms:
            if other.get("element") in ("H", "D"):
                continue
            d = dist(atom["xyz"], other["xyz"])
            if d > 4.5:
                continue
            types = ["general_proximity"]
            if d <= 3.5 and (_hbond(atom, other) or _hbond(other, atom)):
                types.append("hydrogen_bond")
            if d <= 4 and finite(atom.get("formal_charge")) and finite(other.get("formal_charge")) and atom["formal_charge"] * other["formal_charge"] < 0:
                types.append("salt_bridge")
            if atom.get("hydrophobic") is True and other.get("hydrophobic") is True:
                types.append("hydrophobic_contact")
            elements = {atom.get("element"), other.get("element")}
            if d <= 2.8 and elements & METALS and elements & {"N", "O", "S"}:
                types.append("metal_coordination")
            if all(e in VDW for e in elements) and d < .65*(VDW[atom["element"]]+VDW[other["element"]]) and other.get("id") not in atom.get("bonded_atom_ids", []):
                types.append("steric_clash")
            rows.append({"variant": variant, "ligand": ligand, "canonical_residue": canonical,
                         "structure_residue": atom["structure_residue"], "chain": atom["chain"],
                         "amino_acid": atom["amino_acid"], "min_distance": d, "distance_unit": "angstrom",
                         "interaction_types": sorted(types), "protein_atom": atom["name"], "ligand_atom": other["name"],
                         "pose_id": pose_id, "replicate_id": replicate_id, "model_id": model_id, "method": method,
                         "experimentally_reported": annotation.get("experimentally_reported"),
                         "conservation": annotation.get("conservation"), "native_function_importance": annotation.get("native_function_importance"),
                         "protected": annotation.get("protected"), "provenance": provenance})
    return rows


def summarize_contacts(rows, pose_registry):
    """Registry includes ALL sampled poses, including those without contacts.

    Entries: variant, ligand, pose_id, model_id, method, replicate_id. Confidence
    interval is Wilson 95% over independent models, not correlated poses.
    """
    registry = {(r["variant"], r["ligand"], r["pose_id"]): r for r in pose_registry}
    if len(registry) != len(pose_registry):
        raise ValueError("Duplicate pose registry identities")
    grouped = defaultdict(list)
    for row in rows:
        reg = registry.get((row["variant"], row["ligand"], row["pose_id"]))
        if not reg or any(reg.get(k) != row.get(k) for k in ("method", "model_id", "replicate_id")):
            raise ValueError("Contact does not match registered pose/model/method")
        grouped[(row["variant"], row["chain"], row["structure_residue"], row["ligand"])].append(row)
    summary = []
    for key, contacts in sorted(grouped.items()):
        variant, chain, structure_residue, ligand = key
        eligible = [r for r in pose_registry if r["variant"] == variant and r["ligand"] == ligand]
        models = {(r["method"], r["model_id"]) for r in eligible}
        supported_models = {(r["method"], r["model_id"]) for r in contacts}
        p, n = len(supported_models)/len(models), len(models)
        z = 1.959963984540054
        center = (p+z*z/(2*n))/(1+z*z/n)
        half = z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
        by_method = {}
        for method in sorted({r["method"] for r in eligible}):
            universe = {r["model_id"] for r in eligible if r["method"] == method}
            support = {r["model_id"] for r in contacts if r["method"] == method}
            by_method[method] = len(support)/len(universe)
        summary.append({**{k: contacts[0].get(k) for k in ("canonical_residue", "amino_acid", "experimentally_reported", "conservation", "native_function_importance", "protected")},
                        "variant": variant, "chain": chain, "structure_residue": structure_residue, "ligand": ligand,
                        "contact_frequency": len({r["pose_id"] for r in contacts})/len(eligible),
                        "independent_model_frequency": p, "independent_models": n,
                        "contact_interval": [max(0, center-half), min(1, center+half)],
                        "interval_method": "Wilson 95% across independent registered models; independence is supplied by run registry",
                        "method_frequencies": by_method, "methods_supporting": [m for m,f in by_method.items() if f > 0],
                        "min_distance": min(r["min_distance"] for r in contacts), "distance_unit": "angstrom",
                        "interaction_types": sorted({t for r in contacts for t in r["interaction_types"]}),
                        "ligand_atoms": sorted({r["ligand_atom"] for r in contacts}),
                        "protein_atoms": sorted({r["protein_atom"] for r in contacts}),
                        "provenance": [r["provenance"] for r in contacts],
                        "limitations": ["Absence outside the distance cutoff does not establish non-binding", "Chemical typing and contact hypotheses do not measure affinity"]})
    return summary


def classify_residues(residues, *, native_ligands, thresholds):
    """Each residue carries herbicide:{frequency,interval,method_frequencies,n},
    native:{ligand:{same}}, conservation, protected, catalytic, spatial evidence.
    Thresholds are preregistered: reproducibility, native_max, selectivity_min,
    method_tolerance, conservation_max, minimum_models, minimum_methods,
    second_shell_max_angstrom. Missing native evidence never means zero.
    """
    required = ("reproducibility", "native_max", "selectivity_min", "method_tolerance",
                "conservation_max", "minimum_models", "minimum_methods", "second_shell_max_angstrom")
    if not all(finite(thresholds.get(k)) for k in required):
        raise ValueError("All preregistered differential-contact thresholds required")
    if any(not 0 <= thresholds[k] <= 1 for k in ("reproducibility", "native_max", "selectivity_min", "method_tolerance", "conservation_max")):
        raise ValueError("Contact probability thresholds must be in [0, 1]")
    if any(not isinstance(thresholds[k], int) or thresholds[k] < 2 for k in ("minimum_models", "minimum_methods")) or thresholds["second_shell_max_angstrom"] <= 0:
        raise ValueError("At least two models and independent methods and positive spatial cutoff required")
    output = []
    for row in residues:
        result = {**row, "classification": "INSUFFICIENT_EVIDENCE", "nomination_eligible": False, "reasons": []}
        herb = row.get("herbicide", {})
        native = row.get("native", {})
        complete = bool(native_ligands) and all(x in native for x in native_ligands)
        records = [herb] + [native.get(x, {}) for x in native_ligands]
        for r in records:
            interval, methods = r.get("interval"), r.get("method_frequencies", {})
            complete = complete and finite(r.get("frequency")) and 0 <= r.get("frequency", -1) <= 1 and isinstance(interval, list) and len(interval) == 2 and all(finite(x) and 0 <= x <= 1 for x in interval) and interval[0] <= r["frequency"] <= interval[1] and isinstance(r.get("n"), int) and r["n"] >= thresholds["minimum_models"] and len(methods) >= thresholds["minimum_methods"] and all(finite(v) and 0 <= v <= 1 for v in methods.values())
        if row.get("protected") is True or row.get("catalytic") is True:
            result["classification"] = "CATALYTIC_OR_PROTECTED"
        elif not complete or row.get("canonical_residue") is None or row.get("protected") is None or row.get("catalytic") is None or not finite(row.get("conservation")):
            result["reasons"].append("Complete native controls, mapping, protection, conservation and replicate evidence required")
        else:
            max_native = max(native[x]["frequency"] for x in native_ligands)
            result.update({"native_frequencies": {x:native[x]["frequency"] for x in native_ligands},
                           "herbicide_frequency": herb["frequency"], "maximum_native_frequency": max_native,
                           "selectivity_difference": herb["frequency"]-max_native,
                           "herbicide_native_ratio": herb["frequency"]/max_native if max_native > 0 else None})
            unstable = any(max(r["method_frequencies"].values())-min(r["method_frequencies"].values()) > thresholds["method_tolerance"] for r in records)
            strong = herb["frequency"] >= thresholds["reproducibility"] and all(v >= thresholds["reproducibility"] for v in herb["method_frequencies"].values())
            native_strong = max_native >= thresholds["reproducibility"]
            selective = herb["interval"][0]-max(native[x]["interval"][1] for x in native_ligands) >= thresholds["selectivity_min"]
            if unstable:
                result["classification"] = "UNSTABLE_OR_METHOD_DEPENDENT_CONTACT"
            elif strong and native_strong:
                result["classification"] = "SHARED_HERBICIDE_NATIVE_CONTACT"
            elif native_strong or row.get("native_function_importance") == "critical":
                result["classification"] = "NATIVE_CRITICAL_CONTACT"
            elif strong and max_native <= thresholds["native_max"] and selective:
                result["classification"] = "HERBICIDE_SELECTIVE_CONTACT"
            elif not strong and finite(row.get("distance_to_herbicide_contact_angstrom")) and 0 < row["distance_to_herbicide_contact_angstrom"] <= thresholds["second_shell_max_angstrom"] and scientific_provenance(row.get("spatial_evidence")):
                result["classification"] = "SECOND_SHELL_CANDIDATE"
            else:
                result["reasons"].append("Contact selectivity not supported throughout uncertainty interval")
            support = row.get("experimentally_reported") is True or row.get("cross_method_binding_site_support") is True
            result["nomination_eligible"] = result["classification"] in ("HERBICIDE_SELECTIVE_CONTACT", "SECOND_SHELL_CANDIDATE") and row["conservation"] <= thresholds["conservation_max"] and support
            if not result["nomination_eligible"]:
                result["reasons"].append("Protected/native/shared/uncertain contacts, conservation or independent site support prevent nomination")
        output.append(result)
    return output


def validate_mutation_rationale(rationale, residue):
    required = ("mutation", "position_selection", "herbicide_interaction_disrupted", "native_interactions_to_retain",
                "structural_support", "evolutionary_support", "expected_physical_effect", "stability_consequence",
                "uncertainty", "falsifying_experiment")
    errors = [f"missing {x}" for x in required if not rationale.get(x)]
    if not residue.get("nomination_eligible"):
        errors.append("Residue has not passed evidence-driven nomination")
    stability = rationale.get("stability_consequence", {})
    if not isinstance(stability, dict) or not scientific_provenance(stability.get("provenance")):
        errors.append("Stability estimate requires computational or experimental provenance")
    if rationale.get("conservative_substitution") is not True:
        errors.append("Conservative substitution support required")
    return {"eligible": not errors, "errors": errors}
