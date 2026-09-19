"""Prepare, validate, and report a replicated local EPSPS docking/affinity campaign."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
from pathlib import Path

from .evaluate_epsps_campaign import read_campaign
from ..app.validators.mutation_validator import validate_mutation
from ..app.chemistry.contacts import residue_contacts_from_pdb_and_sdf


CONDITIONS = {"native": ("pep", "B"), "herbicide": ("glyphosate", "B"), "native_s3p": ("pep", "C")}
LEGEND = (
    "WT is the unmutated Arabidopsis EPSPS reference (P05466 residues 77-520); mutations use full-sequence numbering. "
    "PEP is phosphoenolpyruvate and S3P is shikimate 3-phosphate, the two native substrates. "
    "Boltz-2 predicts affinity for glyphosate with S3P present, PEP with S3P present, and S3P with PEP present. "
    "Values are predicted pIC50 = 6 - affinity_pred_value, dimensionless; higher means stronger predicted affinity. "
    "Means and brackets summarize two independently seeded runs as mean [minimum, maximum], not confidence intervals. "
    "These model predictions are not measured IC50, Kd, Km, catalytic activity, or validated resistance. "
    "TM-scores, lDDT and RMSDs use all mutant/WT replicate pairs within each of the three prediction contexts; "
    "the table shows worst-case metrics. WT uses an identity reference; separate WT/WT variability controls are reported. "
    "Docking uses cached BioNeMo protein-only structures, four poses per structure, two structures per primary ligand; "
    "S3P is excluded from docking and is not a separate DiffDock ligand in this run. "
    "DiffDock raw confidence is not affinity or probability. New Boltz predictions use query-only MSA and the stored "
    "PubChem chemical states without protonation enumeration. No pIC50 is supplied to the Kd-based retention gate. "
    "Missing Kd uncertainty and folding ddG prevent a complete retention assessment."
)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mutate(sequence, mutation):
    if mutation == "WT":
        return sequence
    _, position, destination = validate_mutation(mutation, sequence, set())
    if not 77 <= position <= 520:
        raise ValueError("Mutation outside modeled mature chain")
    return sequence[:position - 1] + destination + sequence[position:]


def validate_affinity(affinity):
    for required in ("affinity_pred_value", "affinity_probability_binary"):
        if required not in affinity:
            raise ValueError(f"Missing affinity output: {required}")
    for key, value in affinity.items():
        if not key.startswith("affinity_"):
            continue
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
            raise ValueError("Invalid affinity output")
        if key.startswith("affinity_probability") and not 0 <= value <= 1:
            raise ValueError("Invalid binder probability")
    return 6 - affinity["affinity_pred_value"]


def docking_contact_support(poses, mutations):
    result = []
    positions = sorted({int(m[1:-1]) for m in mutations if m != "WT"})
    for position in positions:
        row = {"sequence_position": position}
        for ligand in ("glyphosate", "pep"):
            selected = [pose for pose in poses if pose["mutation"] == "WT" and pose["ligand"] == ligand]
            row[ligand + "_poses"] = len(selected)
            row[ligand + "_contacts"] = sum(position in pose["contacts_uniprot"] for pose in selected)
        result.append(row)
    return result


def prepare(inputs, workflow, contacts, cached, output):
    import numpy as np
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.SeqUtils import seq1
    from rdkit import Chem

    if output.exists():
        raise ValueError("Use a fresh campaign directory")
    previous, groups = read_campaign(cached)
    sequence = "".join(line.strip() for line in (inputs / "P05466.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    if len(sequence) != 520 or previous["sequence"] != sequence:
        raise ValueError("Cached campaign does not match P05466")
    candidates = json.loads((workflow / "candidate_mutations.json").read_text())
    fingerprint = json.loads((workflow / "interaction_fingerprint.json").read_text())
    contact_data = json.loads(contacts.read_text())["records"]
    native = {p for name, row in contact_data.items() for key, values in row.items()
              if key == "s3p_contacts_uniprot" or
              (name.startswith("epsps-native") and key == "primary_ligand_contacts_uniprot") for p in values}
    herbicide = {p for name, row in contact_data.items() if name.startswith("epsps-glyphosate")
                 for p in row["primary_ligand_contacts_uniprot"]}
    protected = set(fingerprint["protected_by_context"])
    allowed = herbicide - native - protected
    if allowed != set(fingerprint["herbicide_selective_mutable"]):
        raise ValueError("Contact fingerprint no longer matches the supplied contact evidence")
    mutations = ["WT"] + [item["mutation"] for item in candidates]
    if len(set(mutations)) != len(mutations) or set(mutations) != {key[0] for key in groups}:
        raise ValueError("Duplicate or mismatched candidate/cached mutation identities")
    for item in candidates:
        _, position, _ = validate_mutation(item["mutation"], sequence, protected)
        if position not in allowed or item["agi"] != "AT2G45300" or item["classification"] != "HERBICIDE_SELECTIVE_MUTABLE":
            raise ValueError("Candidate is not an allowed herbicide-contact mutation")
    ligand_smiles = {name: json.loads((inputs / f"{name}.json").read_text())["PropertyTable"]["Properties"][0]["SMILES"]
                     for name in ("pep", "s3p", "glyphosate")}
    if any(Chem.MolFromSmiles(s) is None for s in ligand_smiles.values()):
        raise ValueError("Invalid ligand SMILES")
    output.mkdir(parents=True)
    (output / "selection").mkdir()
    for filename in ("candidate_mutations.json", "score_packets.json", "interaction_fingerprint.json"):
        shutil.copy2(workflow / filename, output / "selection" / filename)
    shutil.copy2(contacts, output / "contact_report.json")
    records, dock_records = [], []
    base = {**previous, "model": "upstream Boltz-2 2.2.1", "records": [], "conditions": list(CONDITIONS),
            "legend": LEGEND, "source_selection_sha256": digest(workflow / "candidate_mutations.json")}
    for replicate in (1, 2):
        for mutation in mutations:
            full = mutate(sequence, mutation)
            for condition, (primary, binder) in CONDITIONS.items():
                name = f"{mutation}-{condition}-{replicate}"
                request = {"version": 1, "sequences": [
                    {"protein": {"id": "A", "sequence": full[76:], "msa": "empty"}},
                    {"ligand": {"id": "B", "smiles": ligand_smiles[primary]}},
                    {"ligand": {"id": "C", "smiles": ligand_smiles["s3p"]}},
                ], "properties": [{"affinity": {"binder": binder}}]}
                request_path = output / "inputs" / f"replicate-{replicate}" / f"{name}.yaml"
                write_json(request_path, request)
                row = {"record_id": name, "mutation": mutation, "condition": condition,
                       "replicate": replicate, "seed": 42 + replicate,
                       "ligand": "s3p" if binder == "C" else primary,
                       "sequence_sha256": hashlib.sha256(full.encode()).hexdigest(),
                       "request_sha256": digest(request_path)}
                records.append(row)
                base["records"].append({k: row[k] for k in
                                        ("record_id", "mutation", "condition", "replicate", "seed", "sequence_sha256")})
                if condition == "native_s3p":
                    continue
                source = cached / f"{name}.cif"
                structure = MMCIFParser(QUIET=True).get_structure(name, source)
                residues = [r for r in structure[0]["A"] if r.id[0] == " "]
                if "".join(seq1(r.resname) for r in residues) != full[76:]:
                    raise ValueError(f"Cached docking receptor sequence mismatch: {name}")
                if not np.isfinite([atom.coord for r in residues for atom in r]).all():
                    raise ValueError("Non-finite cached receptor coordinates")

                class ProteinOnly(Select):
                    def accept_model(self, model):
                        return model.id == 0

                    def accept_chain(self, chain):
                        return chain.id == "A"

                    def accept_residue(self, residue):
                        return residue.id[0] == " "

                receptor = output / "inputs/receptors" / f"{name}.pdb"
                receptor.parent.mkdir(exist_ok=True)
                writer = PDBIO()
                writer.set_structure(structure)
                writer.save(str(receptor), ProteinOnly())
                dock_records.append({**row, "shard": f"{condition}-{replicate}",
                                     "cached_source": str(source), "source_sha256": digest(source),
                                     "receptor_sha256": digest(receptor)})
    for shard in sorted({row["shard"] for row in dock_records}):
        path = output / "inputs" / f"diffdock-{shard}.csv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["complex_name", "protein_path", "ligand_description", "protein_sequence"])
            writer.writeheader()
            for row in dock_records:
                if row["shard"] == shard:
                    writer.writerow({"complex_name": row["record_id"],
                                     "protein_path": f"/validation/inputs/receptors/{row['record_id']}.pdb",
                                     "ligand_description": ligand_smiles[row["ligand"]], "protein_sequence": ""})
    write_json(output / "campaign/campaign_manifest.json", base)
    manifest = {"target_agi": "AT2G45300", "sequence": sequence, "mutations": mutations,
                "ligands": ligand_smiles, "records": records, "docking_records": dock_records,
                "contact_report_sha256": digest(output / "contact_report.json"), "legend": LEGEND}
    write_json(output / "inputs/input_manifest.json", manifest)
    return {"mutations_including_wt": len(mutations), "boltz_jobs": len(records),
            "docking_complexes": len(dock_records), "expected_docking_poses": 4 * len(dock_records)}


def collect(root, affinity_only=False):
    import numpy as np
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    from rdkit import Chem

    manifest = json.loads((root / "inputs/input_manifest.json").read_text())
    if digest(root / "contact_report.json") != manifest["contact_report_sha256"]:
        raise ValueError("WT contact evidence changed since preparation")
    contact_report = json.loads((root / "contact_report.json").read_text())
    cutoff = contact_report["distance_cutoff_angstrom"]
    affinities, poses = [], []
    for row in manifest["records"]:
        name, replicate = row["record_id"], row["replicate"]
        request_path = root / "inputs" / f"replicate-{replicate}" / f"{name}.yaml"
        if digest(request_path) != row["request_sha256"]:
            raise ValueError("Input changed since preparation")
        request = json.loads(request_path.read_text())
        prediction = root / "boltz" / f"replicate-{replicate}" / f"boltz_results_replicate-{replicate}" / "predictions" / name
        structure_path = prediction / f"{name}_model_0.cif"
        structure = MMCIFParser(QUIET=True).get_structure(name, structure_path)
        observed = "".join(seq1(r.resname) for r in structure[0]["A"] if r.id[0] == " ")
        if observed != mutate(manifest["sequence"], row["mutation"])[76:]:
            raise ValueError(f"Predicted mutation identity mismatch: {name}")
        if not np.isfinite([a.coord for a in structure.get_atoms()]).all():
            raise ValueError("Non-finite Boltz coordinates")
        for entry in request["sequences"][1:]:
            ligand = entry["ligand"]
            expected = Chem.MolFromSmiles(ligand["smiles"]).GetNumHeavyAtoms()
            if sum(a.element not in {"H", "D"} for a in structure[0][ligand["id"]].get_atoms()) != expected:
                raise ValueError("Incomplete ligand atom coverage")
        affinity_path = prediction / f"affinity_{name}.json"
        affinity = json.loads(affinity_path.read_text())
        pic50 = validate_affinity(affinity)
        affinities.append({**row, **affinity, "predicted_pIC50": pic50,
                           "structure": str(structure_path.relative_to(root)), "structure_sha256": digest(structure_path),
                           "affinity_artifact": str(affinity_path.relative_to(root)), "affinity_sha256": digest(affinity_path)})
        shutil.copy2(structure_path, root / "campaign" / f"{name}.cif")
    write_json(root / "affinity_predictions.json", {"records": affinities, "legend": LEGEND})
    campaign_path = root / "campaign/campaign_manifest.json"
    campaign = json.loads(campaign_path.read_text())
    campaign.update(conditions=list(CONDITIONS), legend=LEGEND, records=[
        {k: row[k] for k in ("record_id", "mutation", "condition", "replicate", "seed", "sequence_sha256")}
        for row in manifest["records"]])
    write_json(campaign_path, campaign)
    write_json(root / "affinity_validation.json", {
        "status": "VALIDATED", "count": len(affinities), "validator_sha256": digest(Path(__file__)),
        "legend": "Execution and artifact-integrity validation, not biological accuracy or affinity calibration.",
    })
    if affinity_only:
        return {"validated_affinity_predictions": len(affinities), "docking_validation": "not_requested"}
    for row in manifest["docking_records"]:
        receptor = root / "inputs/receptors" / f"{row['record_id']}.pdb"
        if digest(receptor) != row["receptor_sha256"]:
            raise ValueError("Docking receptor changed since preparation")
        folder = root / "diffdock" / row["shard"] / row["record_id"]
        expected = Chem.MolToSmiles(Chem.MolFromSmiles(manifest["ligands"][row["ligand"]]))
        rows = []
        for path in folder.glob("rank*_confidence*.sdf"):
            match = re.fullmatch(r"rank(\d+)_confidence(-?\d+(?:\.\d+)?)\.sdf", path.name)
            if match is None:
                raise ValueError("Invalid pose filename")
            molecules = list(Chem.SDMolSupplier(str(path), removeHs=True))
            if len(molecules) != 1 or molecules[0] is None or Chem.MolToSmiles(molecules[0]) != expected:
                raise ValueError("Docked ligand identity mismatch")
            xyz = molecules[0].GetConformer().GetPositions()
            if not np.isfinite(xyz).all() or np.ptp(xyz, axis=0).max() <= 0:
                raise ValueError("Invalid docking coordinates")
            rows.append({**row, "rank": int(match[1]), "confidence": float(match[2]),
                         "artifact": str(path.relative_to(root)), "artifact_sha256": digest(path),
                         "contacts_uniprot": residue_contacts_from_pdb_and_sdf(
                             receptor.read_text(), Chem.MolToMolBlock(molecules[0]), cutoff,
                             protein_chain="A", residue_offset=76)})
        rows.sort(key=lambda r: r["rank"])
        if [r["rank"] for r in rows] != [1, 2, 3, 4] or any(a["confidence"] < b["confidence"] for a, b in zip(rows, rows[1:])):
            raise ValueError("Missing, duplicate, or unsorted docking poses")
        poses.extend(rows)
    write_json(root / "docking_predictions.json", {"records": poses, "legend": LEGEND})
    write_json(root / "docking_contact_support.json", {
        "rows": docking_contact_support(poses, manifest["mutations"]), "cutoff_angstrom": cutoff,
        "legend": "WT pose counts with at least one ligand/protein heavy-atom distance within the stated cutoff. "
                  "All ranks are included without confidence filtering. Counts are geometric hypotheses, not binding probabilities. "
                  "PEP is the only native ligand docked; excluded S3P means this is not a complete native-contact protection check.",
    })
    return {"validated_affinity_predictions": len(affinities), "validated_docking_poses": len(poses)}


def report(root, git_commit=None):
    execution_path = root / "execution_manifest.json"
    execution = json.loads(execution_path.read_text()) if execution_path.is_file() else None
    if execution and execution["status"] not in {"COMPUTATION_COMPLETE_ANALYSIS_PENDING", "AVAILABLE_STAGES_COMPLETE"}:
        raise ValueError("Cannot mark an incomplete or failed computation as complete")
    manifest = json.loads((root / "inputs/input_manifest.json").read_text())
    affinities = json.loads((root / "affinity_predictions.json").read_text())["records"]
    docking = json.loads((root / "docking_predictions.json").read_text())["records"]
    assessment = root / "assessment/AT2G45300-glyphosate"
    evidence = json.loads((assessment / "retention_evidence.json").read_text())
    decisions = {r["mutation"]: r["decision"] for r in json.loads((assessment / "function_retention_report.json").read_text())}
    rows = []
    lines = ["# EPSPS Docking, Affinity and Structure Campaign", "",
             "The available EPSPS computational stages completed. Full biological validation remains incomplete.", "",
             "| Mutation | Glyphosate pIC50 | PEP pIC50 | S3P pIC50 | Min qTM / tTM | Min CA lDDT | Max global / site RMSD (A) | Retention decision |",
             "|---|---|---|---|---|---|---|---|"]
    for mutation in manifest["mutations"]:
        row = {"mutation": mutation}
        cells = []
        for ligand in ("glyphosate", "pep", "s3p"):
            data = [a for a in affinities if a["mutation"] == mutation and a["ligand"] == ligand]
            if sorted(a["replicate"] for a in data) != [1, 2]:
                raise ValueError("Incomplete affinity replicate coverage")
            values = [a["predicted_pIC50"] for a in data]
            row.update({f"{ligand}_pic50_mean": statistics.mean(values), f"{ligand}_pic50_min": min(values),
                        f"{ligand}_pic50_max": max(values)})
            cells.append(f"{statistics.mean(values):.3f} [{min(values):.3f}, {max(values):.3f}]")
        metrics = evidence[mutation]["structural_metrics"]
        row.update(metrics)
        row["decision"] = decisions[mutation]
        row["docking_pose_count"] = sum(r["mutation"] == mutation for r in docking)
        rows.append(row)
        def fmt(value):
            return "N/A" if value is None else f"{value:.3f}"
        lines.append("| " + " | ".join([mutation, *cells,
                     f"{fmt(metrics['query_tm_score'])} / {fmt(metrics['target_tm_score'])}", fmt(metrics["alignment_lddt"]),
                     f"{fmt(metrics['global_ca_rmsd_angstrom'])} / {fmt(metrics['active_site_rmsd_angstrom'])}", decisions[mutation]]) + " |")
    lines += ["", "**Table legend.** " + LEGEND, ""]
    contact_path = root / "docking_contact_support.json"
    if contact_path.is_file():
        support = json.loads(contact_path.read_text())
        lines += ["## Independent Docking Contact Support", "",
                  "| Selected WT residue | Glyphosate contact poses / total | PEP contact poses / total |",
                  "|---|---:|---:|"]
        for item in support["rows"]:
            lines.append(f"| {item['sequence_position']} | {item['glyphosate_contacts']} / {item['glyphosate_poses']} | "
                         f"{item['pep_contacts']} / {item['pep_poses']} |")
        lines += ["", f"**Table legend.** Cutoff: {support['cutoff_angstrom']} A. " + support["legend"], ""]
    lines += ["## Unavailable Stages", "",
              "- Other five target systems: validated live complex inputs are not configured.",
              "- GPT-Rosalind evidence agents and LLM judge: live integration/access is not configured; no mock reviews are represented as live.",
              "- Matched Kd prediction intervals, folding ddG and experimental assays: unavailable.",
              "- Two-seed prediction ranges do not establish native-substrate equivalence or absence of herbicide binding.", ""]
    (root / "CAMPAIGN_REPORT.md").write_text("\n".join(lines))
    with (root / "mutation_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(root / "mutation_summary.json", {"rows": rows, "legend": LEGEND})
    (root / "mutation_summary_legend.md").write_text(LEGEND + "\n")
    if execution:
        execution.update(status="AVAILABLE_STAGES_COMPLETE", overall_status="INCOMPLETE_BLOCKED_STAGES",
                         analysis={"git_commit": git_commit, "structural_comparisons": "completed", "retention_gate": decisions,
                                   "combined_table": "mutation_summary.csv"})
        write_json(execution_path, execution)
    return {"decisions": decisions, "affinity_predictions": len(affinities), "docking_poses": len(docking)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("inputs", "workflow", "contacts", "cached", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    col = sub.add_parser("collect")
    col.add_argument("root", type=Path)
    col.add_argument("--affinity-only", action="store_true")
    rep = sub.add_parser("report")
    rep.add_argument("root", type=Path)
    rep.add_argument("--git-commit")
    args = parser.parse_args()
    values = vars(args)
    command = values.pop("command")
    result = prepare(**values) if command == "prepare" else globals()[command](**values)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
