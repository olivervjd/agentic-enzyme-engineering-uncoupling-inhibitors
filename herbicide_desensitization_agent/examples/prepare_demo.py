"""Export only allowlisted, verified scientific artifacts for a read-only local demo."""
from __future__ import annotations

import argparse
import json
import shutil
import statistics
from pathlib import Path

from ..app.backends.calibration import sha, write
from ..app.backends.binding_calibration import build_binding_report
from ..app.registry.loader import TargetRegistry


def inside(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Artifact path escapes the selected run")
    return path


def geometry(path, offset=0, reference=None):
    from Bio.PDB import MMCIFParser, Superimposer
    model = MMCIFParser(QUIET=True).get_structure("model", path)[0]
    if reference:
        ref = MMCIFParser(QUIET=True).get_structure("reference", reference)[0]
        fixed = [r["CA"] for r in ref["A"] if r.id[0] == " " and "CA" in r]
        moving = [r["CA"] for r in model["A"] if r.id[0] == " " and "CA" in r]
        if len(fixed) != len(moving):
            raise ValueError("Viewer overlay requires matching observed residues")
        fitter = Superimposer()
        fitter.set_atoms(fixed, moving)
        fitter.apply(list(model.get_atoms()))
    coords = lambda a: [round(float(x), 3) for x in a.coord]
    protein = [{"position": r.id[1] + offset, "name": r.resname, "xyz": coords(r["CA"])}
               for r in model["A"] if r.id[0] == " " and "CA" in r]
    ligands = [{"chain": c.id, "atoms": [{"element": a.element, "xyz": coords(a)} for a in c.get_atoms()
                                        if a.element not in {"H", "D"}]}
               for c in model if c.id != "A"]
    return {"protein": protein, "ligands": ligands}


def export(calibration, archive, output, model_access=None, binding_dataset=None, workflow_root=None):
    calibration, archive, output = map(Path, (calibration, archive, output))
    report = json.loads((calibration / "calibration_report.json").read_text())
    for item in report["artifacts"]:
        if sha(inside(calibration, item["path"])) != item["sha256"]:
            raise ValueError("Changed calibration artifact: " + item["path"])
    manifest = json.loads((calibration / "calibration_inputs.json").read_text())
    predictions = json.loads((calibration / "predictions.json").read_text())["records"]
    binding = build_binding_report(calibration, binding_dataset)
    structures, summary = {}, []
    for subject in manifest["subjects"]:
        metrics = [r for r in report["comparisons"] if r["subject"] == subject]
        row = {"subject": subject, "mutation": "G96A" if subject == "ecoli_G96A" else "WT",
               "organism": "Arabidopsis" if subject == "arabidopsis_wt" else "E. coli control",
               "tm_min": min(r["query_tm_score"] for r in metrics),
               "rmsd_max": max(r["global_ca_rmsd_angstrom"] for r in metrics), "ligands": {}}
        for ligand in ("glyphosate", "pep", "s3p"):
            values = [r["predicted_pIC50"] for r in predictions if r["subject"] == subject and r["ligand"] == ligand]
            if values:
                row["ligands"][ligand] = {"mean": statistics.mean(values), "min": min(values), "max": max(values), "n": len(values)}
        summary.append(row)
    for prediction in predictions:
        item = {k: prediction[k] for k in ("id", "subject", "context", "seed", "ligand", "predicted_pIC50")}
        offset = manifest["subjects"][prediction["subject"]].get("offset", 0)
        path = inside(calibration, prediction["structure"])
        reference = next(r for r in predictions if r["subject"] == prediction["subject"] and r["context"] == prediction["context"])
        item.update(geometry(path, offset, inside(calibration, reference["structure"])), sha256=sha(path),
                    contacts={k: [p + offset for p in v] for k, v in prediction["contacts"].items()})
        structures[item["id"]] = item
    workflow_directory = Path(workflow_root) if workflow_root else calibration / "workflow"
    run_directory = workflow_directory / "AT2G45300-glyphosate"
    workflow = json.loads((run_directory / "workflow_stages.json").read_text())
    stages = [{"stage": r["stage"], "status": r["status"],
               "reasons": r["detail"].get("reasons", []) if isinstance(r.get("detail"), dict) else
                          [r["detail"]] if isinstance(r.get("detail"), str) else []} for r in workflow]
    diagnostic_review = json.loads((run_directory / "workflow_review.json").read_text()) if (run_directory / "workflow_review.json").exists() else None
    access = json.loads(Path(model_access).read_text()) if model_access else None
    # Only these fields can reach the browser, never raw authentication or environment data.
    public_access = {k: access[k] for k in ("judge", "evidence", "scope", "calibration_review") if k in access} if access else None
    if workflow_root and (workflow_directory / "model_configuration.json").exists():
        configuration = json.loads((workflow_directory / "model_configuration.json").read_text())
        public_access = {role: {k: value[k] for k in ("model", "metadata_accessible", "http_status", "status") if k in value}
                         for role, value in configuration.items()}
        public_access["scope"] = "Recorded diagnostic workflow; candidate approval remains blocked"
    data = {"target": {"id": manifest["target_agi"], "name": "EPSPS", "organism": "Arabidopsis thaliana", "herbicide": "Glyphosate"},
            "targets": [{"id": e.agi, "name": e.protein_name, "herbicide": e.herbicide} for e in TargetRegistry.default().entries],
            "run": {"id": "epsps-calibration-matchedmsa-20260920", "mode": "Recorded run", "predictions": len(predictions),
                    "comparisons": len(report["comparisons"]), "passed": sum(not r["failed_metrics"] for r in report["comparisons"]),
                    "status": "NEEDS_REASSESSMENT", "candidates_advanced": 0},
            "checks": report["checks"], "blocking_reasons": report["blocking_reasons"], "contact_audit": report["contact_audit"],
            "comparisons": report["comparisons"], "thresholds": report["thresholds"], "summary": summary,
            "binding": binding, "workflow": stages, "models": public_access, "workflow_review": diagnostic_review,
            "sources": json.loads((calibration / "evidence.json").read_text())["sources"],
            "archive": json.loads((archive / "mutation_summary.json").read_text()),
            "legends": {"calibration": report["legend"], "structure": (calibration / "calibration_figure_legend.md").read_text(),
                        "archive": (archive / "mutation_summary_legend.md").read_text()},
            "integrity": {"verified_artifacts": len(report["artifacts"]), "calibration_report_sha256": sha(calibration / "calibration_report.json"),
                          "archive_summary_sha256": sha(archive / "mutation_summary.json")}}
    data["run"]["workflow_id"] = workflow_directory.name
    retrieval_path = workflow_directory / "literature/retrieval.json"
    if workflow_root and retrieval_path.exists():
        retrieval = json.loads(retrieval_path.read_text())
        data["literature"] = {k: retrieval[k] for k in ("status", "query", "retrieved_at", "legend")}
        data["sources"].extend({k: s[k] for k in ("id", "url", "scope", "title", "sha256", "evidence_level")}
                               for s in retrieval["sources"])
    output.mkdir(parents=True, exist_ok=True)
    write(output / "results.json", data)
    write(output / "structures.json", structures)
    write(output / "binding_calibration.json", binding)
    shutil.copyfile(calibration / "calibration_structure_checks.png", output / "structure-checks.png")
    print(f"Exported {len(structures)} verified structures; no credentials or runtime secrets included")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model-access", type=Path)
    parser.add_argument("--binding-dataset", type=Path)
    parser.add_argument("--workflow-root", type=Path, help="Show an actual newer diagnostic workflow without altering scientific predictions")
    args = parser.parse_args()
    export(args.calibration, args.archive, args.output, args.model_access, args.binding_dataset, args.workflow_root)


if __name__ == "__main__":
    main()
