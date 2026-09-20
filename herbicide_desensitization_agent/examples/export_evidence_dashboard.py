"""Reanalyse verified recorded runs and export a portable evidence workbench.

Coordinates and recorded tool outputs supply all numerical scientific values.
This command never infers binding, function, or a calibrated interval from a pose.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import math
import shutil
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path

from .prepare_demo import export as export_base, geometry, inside
from ..app.backends.calibration import sha, write

MODEL = "Boltz-2 2.2.1"
CONTACT_LEGEND = ("General proximity is minimum protein–ligand heavy-atom distance <=5.0 Å. "
                  "Frequencies count one vote per independent seed within one ligand context; absence counts as zero. "
                  "A 3/3 frequency is an observed fraction, not a calibrated probability. Chemical hydrogen bonds, "
                  "salt bridges, aromatic and water-mediated contacts require validated atom typing/geometry; "
                  "they are not inferred from proximity. Conservation is unweighted query-identity fraction "
                  "among nongap homolog rows, not a calibrated evolutionary constraint score.")


def read(path):
    return json.loads(Path(path).read_text())


def table_export(root, name, rows, legend):
    write(root / f"{name}.json", {"rows": rows, "legend": legend})
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (root / f"{name}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or ["status"])
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, allow_nan=False) if isinstance(v, (dict, list))
                             else v for k, v in row.items()})
    (root / f"{name}.legend.txt").write_text(legend + "\n")


def contact_records(model, prediction, offset):
    """Retain every contacting heavy atom pair, not only an opaque fingerprint."""
    import numpy as np
    ligand_chains = {"B": "glyphosate" if prediction["context"] == "herbicide" else "pep", "C": "s3p"}
    if prediction["context"] == "s3p_only":
        ligand_chains = {"B": "s3p"}
    rows = []
    for chain, ligand in ligand_chains.items():
        if chain not in model:
            continue
        atoms = [a for a in model[chain].get_atoms() if a.element not in {"H", "D"}]
        for residue in model["A"]:
            if residue.id[0] != " ":
                continue
            pairs = []
            for atom in residue:
                if atom.element in {"H", "D"}:
                    continue
                for other in atoms:
                    distance = float(np.linalg.norm(atom.coord - other.coord))
                    if distance <= 5.0:
                        pairs.append({"protein_atom": atom.name, "ligand_atom": other.name,
                                      "ligand_chain": chain, "distance_angstrom": distance,
                                      "protein_xyz": [float(x) for x in atom.coord],
                                      "ligand_xyz": [float(x) for x in other.coord],
                                      "interaction_type": "general proximity contact"})
            if pairs:
                nearest = min(pairs, key=lambda x: x["distance_angstrom"])
                rows.append({"variant": prediction["subject"], "residue": residue.id[1] + offset,
                             "structure_residue": residue.id[1], "insertion_code": residue.id[2].strip(),
                             "chain": "A", "amino_acid": residue.resname, "ligand": ligand,
                             "pose_id": prediction["id"], "replicate_id": prediction["seed"],
                             "context": prediction["context"], "method": MODEL,
                             "min_distance": nearest["distance_angstrom"],
                             "protein_atom": nearest["protein_atom"], "ligand_atom": nearest["ligand_atom"],
                             "interaction_types": ["general proximity contact"], "atom_pairs": pairs})
    return rows


def conservation(calibration, manifest):
    from Bio import SeqIO
    result = {}
    for subject, info in manifest["subjects"].items():
        rows = ["".join(c for c in str(r.seq) if not c.islower() and c != ".")
                for r in SeqIO.parse(inside(calibration, info["msa"]["path"]), "fasta")]
        for index, aa in enumerate(rows[0]):
            column = [r[index] for r in rows[1:] if r[index] not in "-X"]
            result[(subject, index + 1 + info.get("offset", 0))] = (
                sum(c == aa for c in column) / len(column) if column else None)
    return result


def aggregate_contacts(raw, predictions, manifest, conservation_scores):
    # Only the requested affinity context for each ligand contributes to its matrix.
    desired = {"glyphosate": "herbicide", "pep": "native", "s3p": "native_s3p"}
    by_key = defaultdict(list)
    for row in raw:
        if row["context"] == desired[row["ligand"]]:
            by_key[(row["variant"], row["residue"], row["ligand"])].append(row)
    result = []
    sites = sorted({(r["variant"], r["residue"]) for r in raw if r["context"] in desired.values()})
    for subject, position in sites:
        counts = {ligand: len([p for p in predictions if p["subject"] == subject and p["context"] == context])
                  for ligand, context in desired.items()}
        freq = {ligand: len(by_key[(subject, position, ligand)]) / counts[ligand] if counts[ligand] else None
                for ligand in desired}
        native = [freq[k] for k in ("pep", "s3p") if freq[k] is not None]
        maximum = max(native) if len(native) == 2 else None
        herb = freq["glyphosate"]
        # Current observations are not a preregistered mutation-selection rule.
        classification = "INSUFFICIENT_EVIDENCE"
        if maximum is not None and herb is not None:
            if herb >= 2/3 and maximum >= 2/3:
                classification = "SHARED_HERBICIDE_NATIVE_CONTACT"
            elif maximum >= 2/3:
                classification = "NATIVE_CRITICAL_CONTACT"
            elif 0 < herb < 2/3:
                classification = "UNSTABLE_OR_METHOD_DEPENDENT_CONTACT"
            elif herb >= 2/3:
                classification = "HERBICIDE_SELECTIVE_CONTACT"
        for ligand in desired:
            records = by_key[(subject, position, ligand)]
            minimum = min(records, key=lambda x: x["min_distance"]) if records else None
            result.append({"variant": subject, "residue": position,
                           "structure_residue": position - manifest["subjects"][subject].get("offset", 0),
                           "chain": "A", "ligand": ligand, "contact_frequency": freq[ligand],
                           "min_distance": minimum["min_distance"] if minimum else None,
                           "distance_unit": "angstrom", "interaction_types": ["general proximity contact"] if records else [],
                           "amino_acid": minimum["amino_acid"] if minimum else None,
                           "ligand_atoms": sorted({a["ligand_atom"] for r in records for a in r["atom_pairs"]}),
                           "protein_atoms": sorted({a["protein_atom"] for r in records for a in r["atom_pairs"]}),
                           "methods_supporting": [MODEL] if records else [], "independent_method_count": int(bool(records)),
                           "replicates": counts[ligand], "supporting_replicates": len(records),
                           "fraction_independent_models": freq[ligand], "fraction_poses": freq[ligand],
                           "pose_ids": [r["pose_id"] for r in records], "classification": classification,
                           "herbicide_contact_frequency": herb, "maximum_native_contact_frequency": maximum,
                           "native_contact_frequencies": {k: freq[k] for k in ("pep", "s3p")},
                           "herbicide_selectivity_difference": herb - maximum if herb is not None and maximum is not None else None,
                           "herbicide_native_contact_ratio": herb / maximum if maximum and herb is not None else None,
                           "method_agreement": "Independent matched method unavailable",
                           "contact_uncertainty": "Three-seed observed fraction; interval uncalibrated",
                           "experimental_support": "Not mapped to same-target experimental contacts",
                           "conservation": conservation_scores.get((subject, position)),
                           "protected_status": "Conservatively held: native contact" if maximum and maximum >= 2/3 else "Unreviewed; not eligible",
                           "native_function_importance": "Native proximity is not catalytic evidence",
                           "distance_from_catalytic_residues": None,
                           "known_resistance_evidence": None,
                           "mutation_feasibility": "Not assessed; calibration gate closed"})
    return result


def replicated_geometry(models, predictions, raw):
    """Same-subject/chemical-context atom-name RMSD after a global C-alpha fit."""
    import numpy as np
    from Bio.PDB import Superimposer
    rows = []
    for left, right in itertools.combinations(predictions, 2):
        if (left["subject"], left["context"]) != (right["subject"], right["context"]):
            continue
        first, second = models[left["id"]], models[right["id"]]
        fitter = Superimposer()
        fitter.set_atoms([r["CA"] for r in first["A"] if r.id[0] == " "],
                         [r["CA"] for r in second["A"] if r.id[0] == " "])
        rotation, translation = fitter.rotran
        pocket = {r["structure_residue"] for r in raw if r["pose_id"] in {left["id"], right["id"]}}
        def rmsd(atom_pairs):
            return float(np.sqrt(np.mean([np.sum((a.coord - (b.coord @ rotation + translation))**2)
                                          for a, b in atom_pairs]))) if atom_pairs else None
        sidechains = [(atom, second["A"][residue.id][atom.name]) for residue in first["A"]
                      if residue.id[0] == " " and residue.id[1] in pocket for atom in residue
                      if atom.name not in {"N", "CA", "C", "O", "OXT"} and atom.element not in {"H", "D"}
                      and atom.name in second["A"][residue.id]]
        ligand_values = {}
        for chain in first:
            if chain.id == "A" or chain.id not in second:
                continue
            target = {(a.parent.id, a.name): a for a in second[chain.id].get_atoms()}
            atoms = [a for a in chain.get_atoms() if a.element not in {"H", "D"}]
            pairs = [(a, target[(a.parent.id, a.name)]) for a in atoms if (a.parent.id, a.name) in target]
            ligand_values[chain.id] = rmsd(pairs) if len(pairs) == len(atoms) else None
        contact_sets = [{(r["structure_residue"], r["ligand"]) for r in raw if r["pose_id"] == p["id"]}
                        for p in (left, right)]
        union = contact_sets[0] | contact_sets[1]
        rows.append({"subject": left["subject"], "context": left["context"], "query": left["id"],
                     "reference": right["id"], "kind": "replicate_geometry",
                     "pocket_sidechain_rmsd_angstrom": rmsd(sidechains), "ligand_rmsd_by_chain_angstrom": ligand_values,
                     "contact_jaccard": len(contact_sets[0] & contact_sets[1]) / len(union) if union else None,
                     "pocket_volume_change": None, "assembly_interface_change": None,
                     "method": "Biopython global CA Kabsch fit; same atom names; no symmetry remapping",
                     "limitations": ["Pocket volume and assembly stability not evaluated", "No ligand-symmetry optimization"]})
    return rows


def experimental_geometry(calibration, manifest, models, predictions, raw):
    """Same-sequence homolog contact recovery and pocket atom displacement.

    Author numbering is explicitly matched to the observed-sequence mapping.
    Ligand atom permutations are never guessed from proximity or atom order.
    """
    import numpy as np
    from Bio.PDB import MMCIFParser, Superimposer
    rows = []
    for subject, context, code in (("ecoli_wt", "herbicide", "1G6S"), ("ecoli_G96A", "s3p_only", "1MI4")):
        reference = MMCIFParser(QUIET=True).get_structure(code, calibration / "references" / f"{code}.cif")[0]
        protein = [r for r in reference["A"] if r.id[0] == " " and "CA" in r]
        numbering = manifest["subjects"][subject]["observed_author_residue_numbers"]
        if [r.id[1] for r in protein] != numbering:
            raise ValueError("Experimental author-to-modeled residue mapping differs")
        ligand_residues = {"glyphosate": [r for r in reference.get_residues() if r.resname == "GPJ"],
                           "s3p": [r for r in reference.get_residues() if r.resname == "S3P"]}
        experimental_contacts = {}
        for ligand, residues in ligand_residues.items():
            xyz = np.array([a.coord for r in residues for a in r if a.element not in {"H", "D"}])
            if len(xyz):
                experimental_contacts[ligand] = {i + 1 for i, r in enumerate(protein)
                    if any(np.linalg.norm(xyz - a.coord, axis=1).min() <= 5 for a in r if a.element not in {"H", "D"})}
        pocket = set().union(*experimental_contacts.values())
        for p in predictions:
            if (p["subject"], p["context"]) != (subject, context):
                continue
            predicted = [r for r in models[p["id"]]["A"] if r.id[0] == " " and "CA" in r]
            if len(protein) != len(predicted):
                raise ValueError("Experimental comparison requires matching observed coverage")
            fit = Superimposer()
            fit.set_atoms([r["CA"] for r in protein], [r["CA"] for r in predicted])
            rotation, translation = fit.rotran
            pairs = [(r[a.name], a) for i, (r, q) in enumerate(zip(protein, predicted), 1) if i in pocket for a in q
                     if a.name in r and a.name not in {"N", "CA", "C", "O", "OXT"} and a.element not in {"H", "D"}]
            rmsd = float(np.sqrt(np.mean([np.sum((a.coord - (b.coord @ rotation + translation)) ** 2) for a, b in pairs]))) if pairs else None
            recovery = {}
            for ligand, sites in experimental_contacts.items():
                observed = {r["structure_residue"] for r in raw if r["pose_id"] == p["id"] and r["ligand"] == ligand}
                recovery[ligand] = {"recall": len(sites & observed) / len(sites),
                                    "precision": len(sites & observed) / len(observed) if observed else None,
                                    "experimental_residues": sorted(sites), "predicted_residues": sorted(observed),
                                    "numbering": "Modeled observed sequence; author mapping in chemical context manifest"}
            rows.append({"subject": subject, "context": context, "query": p["id"], "reference": code,
                         "kind": "experimental_pocket_accuracy_homolog", "pocket_sidechain_rmsd_angstrom": rmsd,
                         "predicted_experimental_contact_recovery": recovery, "ligand_rmsd_angstrom": None,
                         "ligand_rmsd_missing_reason": "Chemically validated atom mapping between experimental CCD and generated ligand unavailable; atom-order RMSD prohibited",
                         "method": "Global CA Kabsch transform; matched named pocket side-chain atoms; heavy-atom contacts <=5 Å",
                         "limitations": ["Homolog control is not same-target Arabidopsis validation", "Protonation contexts not experimentally matched"]})
    return rows


def export(calibration, archive, workflow, output):
    started = time.monotonic()
    calibration, archive, workflow, output = map(Path, (calibration, archive, workflow, output))
    data_dir = output / "data"
    results = export_base(calibration, archive, data_dir, workflow_root=workflow)
    manifest = read(calibration / "calibration_inputs.json")
    report = read(calibration / "calibration_report.json")
    predictions = read(calibration / "predictions.json")["records"]
    execution = read(calibration / "execution.json")
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    models, raw = {}, []
    for row in predictions:
        model = parser.get_structure(row["id"], inside(calibration, row["structure"]))[0]
        models[row["id"]] = model
        raw.extend(contact_records(model, row, manifest["subjects"][row["subject"]].get("offset", 0)))
    contacts = aggregate_contacts(raw, predictions, manifest, conservation(calibration, manifest))
    geometries = replicated_geometry(models, predictions, raw) + experimental_geometry(calibration, manifest, models, predictions, raw)
    archive_comparisons_path = archive / "assessment/AT2G45300-glyphosate/structure_comparisons.json"
    archive_comparisons = [{**r, "subject": "archive:" + r["mutation"], "context": r["condition"],
                            "kind": "archive_wt_variability" if r["mutation"] == "WT-control" else "archive_identity_control" if r["mutation"] == "WT" else "archive_mutant_vs_wt",
                            "protocol_scope": "query-only MSA; do not compare to current matched-MSA WT",
                            "source_sha256": sha(archive_comparisons_path)} for r in read(archive_comparisons_path)]
    bindings = []
    for subject in manifest["subjects"]:
        for ligand in ("glyphosate", "pep", "s3p"):
            group = [p for p in predictions if p["subject"] == subject and p["ligand"] == ligand]
            values = [p["predicted_pIC50"] for p in group]
            subset = [r for r in contacts if r["variant"] == subject and r["ligand"] == ligand and r["contact_frequency"]]
            bindings.append({"variant": subject, "ligand": ligand, "role": "herbicide" if ligand == "glyphosate" else "native substrate",
                             "evidence_state": "Insufficient evidence" if group else "Not evaluated",
                             "affinity_metric": "predicted pIC50", "unit": "dimensionless (-log10 IC50/M)",
                             "estimate": statistics.mean(values) if values else None, "interval": None,
                             "observed_range": [min(values), max(values)] if values else None,
                             "interval_description": "No calibrated interval available; observed seed range is descriptive only",
                             "replicates": len(group), "independent_methods": 1 if group else 0,
                             "pose_consistency": [r for r in geometries if r["subject"] == subject and group and r["context"] == group[0]["context"]],
                             "contact_reproducibility": statistics.mean(r["contact_frequency"] for r in subset) if subset else None,
                             "wt_relative_change": None, "method_agreement": "Independent equivalent-context evidence absent",
                             "decision": "INSUFFICIENT_EVIDENCE", "replicate_results": [{"id": p["id"], "seed": p["seed"], "value": p["predicted_pIC50"]} for p in group],
                             "provenance": [p["structure"] for p in group], "limitations": report["blocking_reasons"]})
    # Historical mutant estimates stay separately identified; never compare to current WT protocol.
    historical = read(archive / "mutation_summary.json")["rows"]
    archive_predictions = read(archive / "affinity_predictions.json")["records"]
    viewer_structures = read(data_dir / "structures.json")
    for p in archive_predictions:
        path = inside(archive, p["structure"])
        if sha(path) != p["structure_sha256"]:
            raise ValueError("Archived structure hash changed: " + p["record_id"])
        reference = next(x for x in archive_predictions if x["mutation"] == "WT" and x["condition"] == p["condition"])
        item = {"id": "archive-" + p["record_id"], "subject": "archive:" + p["mutation"],
                "context": p["condition"], "seed": p["seed"], "ligand": p["ligand"],
                "predicted_pIC50": p["predicted_pIC50"], "sha256": p["structure_sha256"],
                "contacts": {}, "wt_reference_id": "archive-" + reference["record_id"],
                "protocol_scope": "Archived query-only MSA, two seeds; not equivalent to current calibration"}
        item.update(geometry(path, 76, inside(archive, reference["structure"])))
        viewer_structures[item["id"]] = item
        target_path = output / "raw/archive" / p["structure"]
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target_path)
    write(data_dir / "structures.json", viewer_structures)
    for row in historical:
        for ligand in ("glyphosate", "pep", "s3p"):
            bindings.append({"variant": "archive:" + row["mutation"], "ligand": ligand,
                             "role": "herbicide" if ligand == "glyphosate" else "native substrate",
                             "evidence_state": "Insufficient evidence", "affinity_metric": "predicted pIC50",
                             "unit": "dimensionless (-log10 IC50/M)", "estimate": row[ligand + "_pic50_mean"],
                             "interval": None, "observed_range": [row[ligand + "_pic50_min"], row[ligand + "_pic50_max"]],
                             "replicates": 2, "pose_consistency": None, "contact_reproducibility": None,
                             "wt_relative_change": row[ligand + "_pic50_mean"] - historical[0][ligand + "_pic50_mean"],
                             "method_agreement": "DiffDock omits S3P: non-equivalent chemical context",
                             "decision": "INSUFFICIENT_EVIDENCE", "legacy_decision": row["decision"],
                             "provenance": "archive/mutation_summary.json", "limitations": ["Query-only MSA; high WT variability; no calibrated interval"]})
    # Preserve true missing evidence rather than back-calibrating the historical candidate set.
    calibration_state = {"status": "INSUFFICIENT_EVIDENCE", "locked": False,
                         "thresholds": report["thresholds"], "threshold_status": "Historical screening cutoffs, not calibrated acceptance thresholds",
                         "reasons": report["blocking_reasons"], "candidate_design_enabled": False}
    from ..app.evidence_system import evaluate_candidate
    decisions = []
    for row in historical[1:]:
        assessment = evaluate_candidate({"variant": row["mutation"]}, None)
        assessment.update(mutation=row["mutation"], legacy_decision=row["decision"],
                          legacy_protocol="query-only MSA, two seeds; not equivalent to current WT")
        decisions.append(assessment)
    mutations = [{"mutation": r["mutation"], "classification": "INSUFFICIENT_EVIDENCE",
                  "herbicide_contacts_disrupted": None, "native_contacts_affected": None,
                  "differential_contact_evidence": "Historical hypothesis; current WT calibration does not support nomination",
                  "conservation": None, "ddg": None, "ddg_unit": "kcal/mol",
                  "structure_retention": {k: v for k, v in r.items() if k.endswith("score") or k.endswith("angstrom") or k.startswith("alignment_")},
                  "binding_evidence": "Uncalibrated predicted pIC50; no direct Kd",
                  "uncertainty": "Cannot isolate mutation effect from WT model variation",
                  "recommendation": "INSUFFICIENT_EVIDENCE", "legacy_decision": r["decision"],
                  "rationale": {"why_selected": "Archived contact-only hypothesis; not renominated",
                                "interaction_to_disrupt": "Unvalidated herbicide proximity",
                                "native_interactions_to_retain": ["PEP with S3P", "S3P with PEP"],
                                "expected_effect": "Not established", "stability_consequence": None,
                                "falsification": "Matched biochemical assay finds native activity loss or no inhibitor weakening"}}
                 for r in historical[1:]]
    controls = [{"control": "Arabidopsis WT repeatability", "expected_outcome": "Stable repeated structures",
                 "predicted_outcome": report["checks"]["wt_repeatability"], "replicate_variability": [r for r in report["comparisons"] if r["kind"] == "repeatability" and r["subject"] == "arabidopsis_wt"],
                 "experimental_agreement": "Same-target experimental structure unavailable", "pass_fail": "PASS_WITHIN_RECORDED_STRUCTURE_SCOPE",
                 "implication": "Does not calibrate affinity or establish enzyme function"}]
    for row in report["control_directions"]:
        controls.append({"control": "E. coli G96A " + row["ligand"], "expected_outcome": "Literature reports weakened inhibitor/substrate response; kinetic endpoints",
                         "predicted_outcome": row, "replicate_variability": row["seed_values"],
                         "experimental_agreement": "NON_EQUIVALENT_ENDPOINTS", "pass_fail": "INSUFFICIENT_EVIDENCE",
                         "implication": "Directional diagnostic only; cannot calibrate pIC50 using Km or Ki"})
    for name in ("Same-target resistant control", "Same-target catalytic-loss control", "Same-target neutral control", "Apo control", "Redocking and cross-docking"):
        controls.append({"control": name, "expected_outcome": "Predeclared validated outcome required", "predicted_outcome": None,
                         "replicate_variability": None, "experimental_agreement": None, "pass_fail": "NOT_EVALUATED",
                         "implication": "Required baseline evidence remains missing"})
    contexts = []
    for p in predictions:
        info = manifest["subjects"][p["subject"]]
        contexts.append({"id": p["id"], "protein_sequence": info["sequence"], "isoform": "Recorded input sequence; isoform not independently curated",
                         "variant": p["subject"], "modeled_residue_range": [1 + info.get("offset", 0), len(info["sequence"]) + info.get("offset", 0)],
                         "chain_mapping": {"A": "protein", "B": "s3p" if p["context"] == "s3p_only" else ("glyphosate" if p["context"] == "herbicide" else "pep"),
                                           "C": None if p["context"] == "s3p_only" else "s3p"},
                         "residue_offset": info.get("offset", 0), "oligomeric_state": "Modeled monomer",
                         "subcellular_context": "plastid" if p["subject"] == "arabidopsis_wt" else "E. coli homolog control",
                         "chemical_state": manifest["chemical_audit"], "ph": None, "ionic_strength": None,
                         "metals": "Not explicitly represented", "catalytic_waters": "Not represented", "cofactors": "Not curated",
                         "structure_model": MODEL, "affinity_model": MODEL, "msa": info["msa"], "seed": p["seed"],
                         "sampling": {"recycling_steps": 3, "sampling_steps": 100, "diffusion_samples": 1, "sampling_steps_affinity": 200, "diffusion_samples_affinity": 3},
                         "provenance": {"request_sha256": p["request_sha256"], "structure_sha256": sha(inside(calibration, p["structure"]))}})
    agent_specs = [
        ("Registry", "Define sequence, ligands, native function and context", "COMPLETED", "Target registry 1.1", ["Canonical registry", "calibration_inputs.json"]),
        ("Literature evidence", "Curate experimental structures and controls", "AWAITING_EVIDENCE", "RCSB primary structures; curated citations", ["1G6S", "1MI4"]),
        ("Chemical state", "Validate ligand identity, stereo, charge and protonation", "AWAITING_EVIDENCE", "RDKit", ["chemical_state_audit.json"]),
        ("Structure", "Replicated WT and homolog-control structures", "COMPLETED", MODEL, ["24 recorded structures; three seeds"]),
        ("Docking", "Independent replicated ligand poses and clusters", "NON_EQUIVALENT_CONTEXT", "DiffDock-L", ["Historical protein-only poses"]),
        ("Affinity", "Ligand-specific estimates and uncertainty", "AWAITING_EVIDENCE", MODEL, ["Predicted pIC50; no calibrated intervals"]),
        ("Contact mapping", "Residue and atom contact maps", "COMPLETED", "Biopython / NumPy heavy-atom distances", ["residue_interactions.json", "pose_contacts.json"]),
        ("Calibration", "WT variability and control response", "AWAITING_EVIDENCE", "TM-align / CA lDDT / control checks", ["model_calibration.json"]),
        ("Mutation design", "Conservative single-substitution hypotheses", "BLOCKED", "Deterministic eligibility gates", ["No new nominations"]),
        ("Stability", "Folding delta-delta-G and assembly integrity", "NOT_EVALUATED", "Not configured", []),
        ("Function retention", "Native-substrate and herbicide gates", "AWAITING_EVIDENCE", "Deterministic evidence engine", ["decisions.json"]),
        ("Evidence review", "Check claims against tool and experimental outputs", "COMPLETED", "Deterministic provenance review", ["limitations.json"]),
        ("Dashboard", "Display validated evidence without changing decisions", "COMPLETED", "Vite / Three.js", ["Interactive dashboard; static report; figures"]),
        ("Learning loop", "Ingest measured assay results without automatic approval", "AWAITING_EVIDENCE", "Governed assay ingestion", ["Experimental-validation plan"]),
    ]
    agents = [{"agent": name, "purpose": purpose, "inputs": inputs, "tools_models": method, "outputs": outputs,
               "status": status, "warnings": [] if status == "COMPLETED" else ["Incomplete evidence; see limitations"],
               "runtime": None, "runtime_unit": "seconds", "provenance": "Recorded scientific run + current deterministic analysis"}
              for name, purpose, status, method, outputs in agent_specs for inputs in [["Verified raw inputs and prior validated stage"]]]
    limitations = report["blocking_reasons"] + [
        "No direct Kd, folding delta-delta-G, pocket-volume or interface-stability measurements.",
        "No new mutation satisfies the calibrated screen. Archived threshold failures do not isolate a mutation effect.",
        "Apo, neutral, catalytic-loss, decoy, same-target resistant controls and matched-context independent docking are incomplete.",
        "No candidate MD, free-energy or biochemical validation has been performed.",
        "Distance-based contacts do not establish hydrogen bonding, catalysis, binding or resistance.",
        "No legal eligibility is asserted: nucleotide edits, crop, jurisdiction and applicable framework require separate review.",
        "Other registered protein–herbicide pairs are not evaluated; no EPSPS conclusions transfer to them."]
    plan = [
        {"step": 1, "action": "Curate WT and same-target resistant, catalytic-loss and neutral controls with primary sources; confirm sequence and numbering.", "readout": "Context-complete baseline", "gate": "No mutation design before baseline acceptance"},
        {"step": 2, "action": "Fix assay pH, ionic conditions, ligand protonation/stereo, substrates and biological assembly; pre-register seeds and thresholds before candidates.", "readout": "Immutable protocol and calibration digest", "gate": "No post hoc threshold changes"},
        {"step": 3, "action": "Measure matched WT/control inhibitor response and native substrate kinetics with biological replicates and reported uncertainty.", "readout": "Ki or IC50 under specified substrate conditions; Km and kcat separately", "gate": "Do not relabel kinetic endpoints as Kd"},
        {"step": 4, "action": "Run matched S3P-aware independent docking, redocking/cross-docking and apo/substrate controls; calibrate contact and pose reproducibility.", "readout": "Ligand-specific poses, clusters, control recovery", "gate": "Independent method agreement required"},
        {"step": 5, "action": "Only after calibration, test minimal single substitutions supported by differential contacts and conservation, then fold/assembly stability and free-energy/MD for shortlisted candidates.", "readout": "Matched mutant-WT effects with complete uncertainty", "gate": "Every native substrate must remain within the preregistered equivalence interval"},
        {"step": 6, "action": "Biochemically validate shortlisted mutants and subsequently confirm phenotype in an appropriate authorised plant experiment.", "readout": "Retained native activity and reduced inhibitor effect", "gate": "Experimentally validated only with traceable measured evidence"},
    ]
    from rdkit import Chem
    from rdkit.Chem import Draw
    ligands = []
    for name, smiles in manifest["ligands"].items():
        mol = Chem.MolFromSmiles(smiles)
        drawer = Draw.MolDraw2DSVG(340, 190)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        ligands.append({"id": name, "name": name, "canonical_identifier": Chem.MolToInchiKey(mol),
                        "smiles": Chem.MolToSmiles(mol), "svg": drawer.GetDrawingText(),
                        "chemical_state": manifest["chemical_audit"][name], "experimental_references": results["sources"]})
    tables = {"binding": bindings, "residue_interactions": contacts, "mutation_selection": mutations,
              "model_calibration": controls, "agent_execution": agents}
    payload = {"schema_version": "1.0", "objective": "Minimal point mutations that reduce herbicide interference while preserving native function",
               "tables": tables, "decisions": decisions, "context_records": contexts, "limitations": limitations,
               "validation_plan": plan, "ligands": ligands, "calibration": calibration_state,
               "structural_comparisons": report["comparisons"] + geometries + archive_comparisons,
               "contact_definitions": CONTACT_LEGEND, "legends": {"contacts": CONTACT_LEGEND,
                   "binding": "Predicted pIC50 is dimensionless; no conversion to Kd/Ki/DeltaG. Ranges are not uncertainty intervals.",
                   "structure": "WT-vs-WT variation is displayed with homolog accuracy; historical candidate runs use a different protocol."},
               "conclusion": {"directly_measured": "Published E. coli crystallographic references; no new assays",
                   "computationally_predicted": "24 recorded Boltz complexes, ligand-specific pIC50 and recalculated coordinate contacts",
                   "methods_agree_on": "WT structural repeatability; independent equivalent-context binding agreement unavailable",
                   "uncertain": limitations, "nominations": "No new nominations; calibration incomplete", "next_experiment": plan[2]},
               "runtime_seconds": time.monotonic() - started}
    write(data_dir / "evidence_system.json", payload)
    table_legends = {"binding": payload["legends"]["binding"], "residue_interactions": CONTACT_LEGEND,
                     "mutation_selection": "Archived hypotheses, not new recommendations. RMSD in angstrom; TM-score/lDDT/coverage dimensionless; folding delta-delta-G in kcal/mol is missing. A legacy failure does not establish a calibrated mutant effect.",
                     "model_calibration": "Controls retain original species, chemical context and physical endpoint. Seed ranges are descriptive, not confidence intervals. Missing control outcomes block calibration.",
                     "agent_execution": "Runtime is seconds where recorded; null means unrecorded. Completed computation is not completed scientific validation. Model commentary cannot replace numerical tools."}
    payload["table_legends"] = table_legends
    write(data_dir / "evidence_system.json", payload)
    for name, rows in tables.items():
        table_export(data_dir, name, rows, table_legends[name])
    for name, value in (("pose_contacts", raw), ("decisions", decisions), ("chemical_contexts", contexts),
                        ("structural_comparisons", payload["structural_comparisons"]), ("limitations", limitations), ("experimental_validation_plan", plan)):
        write(data_dir / f"{name}.json", value)
    (output / "limitations.md").write_text("# Limitations\n\n" + "\n".join("- " + x for x in limitations) + "\n")
    (output / "experimental-validation-plan.md").write_text("# Experimental validation plan\n\n" + "\n\n".join(f"{r['step']}. {r['action']} Readout: {r['readout']} Gate: {r['gate']}" for r in plan))
    # Copy only evidence artifacts; never runtime credential files or whole home directories.
    for item in report["artifacts"]:
        source = inside(calibration, item["path"])
        target = output / "raw" / "calibration" / item["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    for name in ("calibration_inputs.json", "calibration_report.json", "predictions.json", "evidence.json", "execution.json"):
        shutil.copyfile(calibration / name, output / "raw/calibration" / name)
    target = output / "raw/archive"
    target.mkdir(parents=True, exist_ok=True)
    for name in ("mutation_summary.json", "affinity_predictions.json", "docking_predictions.json", "execution_manifest.json"):
        shutil.copyfile(archive / name, target / name)
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "biopython", "rdkit", "tmtools", "matplotlib")}
    write(output / "workflow-manifest.json", {"schema_version": "1.0", "status": "INSUFFICIENT_EVIDENCE", "scientific_claim": "No resistance claim",
          "nodes": agents, "decisions": decisions, "calibration": calibration_state, "analysis_tool_versions": versions,
          "source_calibration_sha256": sha(calibration / "calibration_report.json"), "source_archive_sha256": sha(archive / "mutation_summary.json"),
          "artifacts": [{"path": str(p.relative_to(output)), "sha256": sha(p)} for p in sorted(output.rglob("*")) if p.is_file()]})
    print(f"Exported {len(bindings)} binding rows, {len(contacts)} residue-ligand rows, {len(raw)} pose contacts; no new mutations nominated")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("calibration", "archive", "workflow", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    export(args.calibration, args.archive, args.workflow, args.output)


if __name__ == "__main__":
    main()
