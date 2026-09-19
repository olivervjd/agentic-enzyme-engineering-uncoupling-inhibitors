"""Validate real DiffDock-L poses and upstream Boltz-2 affinity artifacts."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


def validate(root: Path) -> dict:
    import numpy as np
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    from rdkit import Chem

    inputs = json.loads((root / "inputs/input_manifest.json").read_text())
    summary = {"target": "AT2G45300", "mutation": "WT", "tools": {}, "ligands": {}}
    summary["diffdock_runtime"] = json.loads((root / "diffdock_runtime.json").read_text())
    for name in ("pep", "glyphosate"):
        expected = Chem.MolToSmiles(Chem.MolFromSmiles(inputs["ligands"][name]))
        poses = []
        for path in sorted((root / "diffdock" / f"epsps-{name}").glob("rank*_confidence*.sdf")):
            match = re.fullmatch(r"rank(\d+)_confidence(-?\d+(?:\.\d+)?)\.sdf", path.name)
            if match is None:
                raise ValueError(f"Unrecognized DiffDock output filename: {path.name}")
            molecule = Chem.SDMolSupplier(str(path), removeHs=True)[0]
            if molecule is None or Chem.MolToSmiles(molecule) != expected:
                raise ValueError(f"Docked ligand chemistry mismatch: {path.name}")
            xyz = molecule.GetConformer().GetPositions()
            if not np.isfinite(xyz).all() or np.ptp(xyz, axis=0).max() <= 0:
                raise ValueError(f"Invalid pose coordinates: {path.name}")
            confidence = float(match[2])
            if not math.isfinite(confidence):
                raise ValueError("Non-finite DiffDock score")
            poses.append({"rank": int(match[1]), "confidence": confidence,
                          "heavy_atoms": molecule.GetNumHeavyAtoms(), "artifact": str(path.relative_to(root))})
        if sorted(pose["rank"] for pose in poses) != [1, 2, 3, 4]:
            raise ValueError(f"Missing or duplicate DiffDock ranks for {name}")
        poses.sort(key=lambda pose: pose["rank"])
        if any(a["confidence"] < b["confidence"] for a, b in zip(poses, poses[1:])):
            raise ValueError("DiffDock ranks are not ordered by confidence")

        predictions = root / "boltz/boltz_results_boltz/predictions" / f"epsps-{name}"
        affinity_path = predictions / f"affinity_epsps-{name}.json"
        affinity = json.loads(affinity_path.read_text())
        value = float(affinity["affinity_pred_value"])
        probability = float(affinity["affinity_probability_binary"])
        if not math.isfinite(value) or not math.isfinite(probability) or not 0 <= probability <= 1:
            raise ValueError("Invalid Boltz-2 affinity output")
        structure_path = predictions / f"epsps-{name}_model_0.cif"
        structure = MMCIFParser(QUIET=True).get_structure(name, structure_path)
        residues = [r for r in structure[0]["A"] if r.id[0] == " "]
        request = json.loads((root / "inputs/boltz" / f"epsps-{name}.yaml").read_text())
        sequence = "".join(seq1(residue.resname) for residue in residues)
        if sequence != request["sequences"][0]["protein"]["sequence"]:
            raise ValueError("Boltz-2 output protein sequence mismatch")
        if not np.isfinite(np.array([a.coord for a in structure.get_atoms()])).all():
            raise ValueError("Non-finite Boltz-2 coordinates")
        for chain, ligand_name in (("B", name), ("C", "s3p")):
            expected_atoms = Chem.MolFromSmiles(inputs["ligands"][ligand_name]).GetNumHeavyAtoms()
            observed = sum(a.element not in {"H", "D"} for a in structure[0][chain].get_atoms())
            if expected_atoms != observed:
                raise ValueError(f"Boltz-2 ligand atom coverage mismatch for {chain}")
        confidence = json.loads((predictions / f"confidence_epsps-{name}_model_0.json").read_text())
        summary["ligands"][name] = {
            "diffdock_poses": poses, "boltz_affinity": affinity, "predicted_pIC50": 6 - value,
            "boltz_confidence": confidence, "protein_residues": len(residues),
            "structure": str(structure_path.relative_to(root)),
            "affinity_artifact": str(affinity_path.relative_to(root)),
        }
    summary["tools"] = {"diffdock_l": "execution-and-output-validation-passed",
                         "boltz2_structure_and_affinity": "execution-and-output-validation-passed"}
    summary["legend"] = (
        inputs["legend"] + " DiffDock confidence is the raw pose score rounded to two decimals in filenames; "
        "higher is better within a run, but this is not affinity or a calibrated probability. "
        "Boltz affinity_pred_value is predicted log10(IC50 in micromolar), lower means stronger predicted binding; "
        "predicted pIC50 = 6 - affinity_pred_value is dimensionless. Neither is a measured Kd. "
        "affinity_probability_binary is a model-predicted binder probability, not an experimental success rate. "
        "Structure confidence is not a structural similarity measurement. "
        "Artifact validation checks data integrity, not pose accuracy, steric validity, or affinity calibration. "
        "One seed and WT only are tested; there are no mutation effects, confidence intervals, "
        "substrate-equivalence conclusions, or claims of experimentally verified affinity."
    )
    (root / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# Local Docking and Affinity Smoke Test", "", "Both tools passed execution and artifact validation.", "",
             "DiffDock-L device: `" + summary["diffdock_runtime"]["device"] + "`.", "",
             "| Protein | Ligand | Docked poses | Top pose score | Predicted log10(IC50/uM) | Predicted pIC50 | Binder probability |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for name, result in summary["ligands"].items():
        lines.append(f"| EPSPS WT | {name} | {len(result['diffdock_poses'])} | "
                     f"{result['diffdock_poses'][0]['confidence']:.2f} | "
                     f"{result['boltz_affinity']['affinity_pred_value']:.4f} | {result['predicted_pIC50']:.4f} | "
                     f"{result['boltz_affinity']['affinity_probability_binary']:.4f} |")
    lines += ["", "**Table legend.** " + summary["legend"], "",
              "These WT predictions are not yet a calibrated affinity baseline for mutation selection.", ""]
    (root / "VALIDATION_REPORT.md").write_text("\n".join(lines))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    print(json.dumps(validate(parser.parse_args().root)["tools"], indent=2))
