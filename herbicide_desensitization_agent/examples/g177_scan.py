"""Explicit G177 saturation scan; does not alter automatic nomination rules."""
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import statistics
import subprocess
import time

AMINO_ACIDS = dict(zip("ACDEFGHIKLMNPQRSTVWY", (
    "Alanine", "Cysteine", "Aspartate", "Glutamate", "Phenylalanine", "Glycine",
    "Histidine", "Isoleucine", "Lysine", "Leucine", "Methionine", "Asparagine",
    "Proline", "Glutamine", "Arginine", "Serine", "Threonine", "Valine",
    "Tryptophan", "Tyrosine")))
SEEDS = (211, 223, 227)
STATES = {
    "pep": "C=C(OP(=O)([O-])[O-])C(=O)[O-]",
    "s3p": "O=C([O-])C1=C[C@@H](OP(=O)([O-])[O-])[C@@H](O)[C@H](O)C1",
    "glyphosate": "O=C([O-])C[NH2+]CP(=O)([O-])[O-]",
}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def mutate(sequence, aa):
    if len(sequence) != 520 or sequence[176] != "G" or aa not in AMINO_ACIDS:
        raise ValueError("Expected P05466 canonical G177 and a standard amino acid")
    mutant = sequence[:176] + aa + sequence[177:]
    assert sum(a != b for a, b in zip(sequence, mutant)) == (aa != "G")
    assert mutant[76:][100] == aa
    return mutant


def replace_query(a3m, sequence, mutant):
    lines = a3m.splitlines(keepends=True)
    headers = [i for i, line in enumerate(lines) if line.startswith(">")]
    if len(headers) < 2:
        raise ValueError("A homolog alignment, not a query-only MSA, is required")
    start, end = headers[:2]
    if "".join(x.strip() for x in lines[start + 1:end]) != sequence:
        raise ValueError("MSA query does not equal mature WT sequence")
    return "".join(lines[:start + 1]) + mutant + "\n" + "".join(lines[end:])


def prepare(root, donor):
    from boltz.data.parse.schema import standardize
    from rdkit import Chem
    if importlib.metadata.version("boltz") != "2.2.1":
        raise ValueError("This campaign pins Boltz 2.2.1")
    if root.exists():
        raise FileExistsError("Use a fresh campaign directory")
    source = json.loads((donor / "calibration_inputs.json").read_text())
    full = source["sequence"]
    if hashlib.sha256(full.encode()).hexdigest() != "ac2d95791932eca90a6986cbe18b4d9ab948b08d42fb68b871b9b777d8b51bcf":
        raise ValueError("Canonical sequence differs from repository baseline")
    subject = source["subjects"]["arabidopsis_wt"]
    assert subject["sequence"] == full[76:]
    msa_path = donor / subject["msa"]["path"]
    if sha(msa_path) != subject["msa"]["sha256"]:
        raise ValueError("Donor alignment integrity mismatch")
    a3m = msa_path.read_text()
    manifest = {"created_at": now(), "target": "Arabidopsis EPSPS P05466 / AT2G45300",
        "canonical_position": 177, "mature_position": 101, "mature_range": [77, 520],
        "wild_type": "G", "full_sequence": full, "donor_manifest_sha256": sha(donor / "calibration_inputs.json"),
        "source_msa": subject["msa"], "seeds": list(SEEDS), "records": [],
        "protocol": {"boltz": "2.2.1", "recycling_steps": 3, "sampling_steps": 100,
            "diffusion_samples": 1, "sampling_steps_affinity": 200, "diffusion_samples_affinity": 3},
        "scope": "User-requested saturation scan at protected G177, not automatic mutation nomination",
        "msa_policy": "Identical WT homolog rows for every variant; replace only query with exact mutant sequence",
        "limitations": ["Exploratory model predictions, not measured Kd or proof of retained catalytic activity",
            "Affinity parser neutralizes scored ligands; requested and encoded chemistry are both recorded",
            "PEP is a substrate: model IC50-scale outputs are not measured PEP affinities or Km",
            "Protein-mutation effects are not calibrated by this scan",
            "Three seeds measure sampling variation, not calibrated confidence intervals",
            "Ligand-bound cofolding, not a sequential apo-fold and independent binding experiment",
            "No new external MSA search; unchanged homolog background isolates query substitution"],
        "chemical_states": {}}
    root.mkdir(parents=True)
    (root / "msa").mkdir()
    for ligand in ("pep", "glyphosate"):
        encoded = standardize(STATES[ligand])
        manifest["chemical_states"][ligand] = {
            "requested_smiles": STATES[ligand], "encoded_smiles": encoded,
            "requested_charge": Chem.GetFormalCharge(Chem.MolFromSmiles(STATES[ligand])),
            "encoded_charge": Chem.GetFormalCharge(Chem.MolFromSmiles(encoded)),
            "s3p_smiles": STATES["s3p"], "s3p_charge": -3}
    fasta = []
    for aa in "G" + "".join(x for x in AMINO_ACIDS if x != "G"):
        name = "WT" if aa == "G" else f"G177{aa}"
        mature = mutate(full, aa)[76:]
        fasta.append(f">{name}|P05466_77-520|canonical177={aa}|model101={aa}\n{mature}\n")
        msa = root / "msa" / f"{name}.a3m"
        msa.write_text(replace_query(a3m, full[76:], mature))
        for seed in SEEDS:
            for ligand in ("pep", "glyphosate"):
                rid = f"{name}-{ligand}-{seed}"
                request = {"version": 1, "sequences": [
                    {"protein": {"id": "A", "sequence": mature, "msa": str(msa)}},
                    {"ligand": {"id": "B", "smiles": STATES[ligand]}},
                    {"ligand": {"id": "C", "smiles": STATES["s3p"]}}],
                    "properties": [{"affinity": {"binder": "B"}}]}
                path = root / "inputs" / f"seed-{seed}" / f"{rid}.yaml"
                save(path, request)
                manifest["records"].append({"id": rid, "variant": name, "amino_acid": aa,
                    "ligand": ligand, "seed": seed, "request": str(path.relative_to(root)),
                    "request_sha256": sha(path), "msa": str(msa.relative_to(root)), "msa_sha256": sha(msa),
                    "mature_sequence": mature})
    assert len(manifest["records"]) == 120
    (root / "variants.fasta").write_text("".join(fasta))
    save(root / "manifest.json", manifest)
    save(root / "status.json", {"state": "PREPARED", "expected": 120, "completed": 0, "updated_at": now()})
    print("Prepared 19 substitutions + WT, 120 cofolding/affinity jobs", flush=True)


def collect(root):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    manifest = json.loads((root / "manifest.json").read_text())
    results, errors = [], []
    for row in manifest["records"]:
        folder = root / "boltz" / f"seed-{row['seed']}" / f"boltz_results_seed-{row['seed']}" / "predictions" / row["id"]
        affinity = folder / f"affinity_{row['id']}.json"
        cif = folder / f"{row['id']}_model_0.cif"
        confidence = folder / f"confidence_{row['id']}_model_0.json"
        if not all(p.is_file() for p in (affinity, cif, confidence)):
            continue
        try:
            raw = json.loads(affinity.read_text())
            conf = json.loads(confidence.read_text())
            value = float(raw["affinity_pred_value"])
            probability = float(raw["affinity_probability_binary"])
            if not math.isfinite(value) or not math.isfinite(probability) or not 0 <= probability <= 1:
                raise ValueError("Invalid affinity value")
            structure = MMCIFParser(QUIET=True).get_structure(row["id"], str(cif))
            sequence = "".join(seq1(r.resname) for r in structure[0]["A"] if r.id[0] == " ")
            if sequence != row["mature_sequence"] or set(c.id for c in structure[0]) != {"A", "B", "C"}:
                raise ValueError("Output sequence or ligand chain mismatch")
            results.append({**{k: row[k] for k in ("id", "variant", "amino_acid", "ligand", "seed")},
                "affinity_pred_value": value, "affinity_probability_binary": probability,
                "complex_plddt": conf.get("complex_plddt"), "ligand_iptm": conf.get("ligand_iptm"),
                "confidence_score": conf.get("confidence_score"), "raw_affinity": raw,
                "artifacts": [{"path": str(p.relative_to(root)), "sha256": sha(p)} for p in (affinity, cif, confidence)]})
        except (ValueError, KeyError, OSError) as exc:
            errors.append({"id": row["id"], "error": str(exc)})
    rows = []
    for aa in "G" + "".join(x for x in AMINO_ACIDS if x != "G"):
        row = {"amino_acid": aa, "name": AMINO_ACIDS[aa], "variant": "WT" if aa == "G" else f"G177{aa}"}
        for ligand in ("pep", "glyphosate"):
            matches = [r for r in results if r["amino_acid"] == aa and r["ligand"] == ligand]
            values = [r["affinity_pred_value"] for r in matches]
            row[ligand] = {"n": len(values), "mean": statistics.mean(values) if values else None,
                "min": min(values) if values else None, "max": max(values) if values else None,
                "sd": statistics.stdev(values) if len(values) > 1 else None,
                "mean_binding_probability": statistics.mean(r["affinity_probability_binary"] for r in matches) if matches else None}
        rows.append(row)
    for row in rows:
        for ligand in ("pep", "glyphosate"):
            row[ligand]["delta_vs_wt"] = (row[ligand]["mean"] - rows[0][ligand]["mean"]
                if row[ligand]["n"] == rows[0][ligand]["n"] == 3 else None)
    report = {"updated_at": now(), "complete": len(results) == 120 and not errors,
        "completed": len(results), "expected": 120, "errors": errors, "rows": rows, "replicates": results,
        "metric": "Boltz affinity_pred_value; log10(IC50 / micromolar); lower predicts stronger binding",
        "limitations": manifest["limitations"]}
    save(root / "results.json", report)
    lines = ["# G177 saturation scan — Boltz2 predictions", "",
        f"Validated results: {len(results)}/120. All complexes include S3P. Canonical residue 177 = mature residue 101.", "",
        "Scores are Boltz affinity_pred_value, log10(IC50 / µM), reported as mean ± sample SD over three seeds; lower predicts stronger binding. These are uncalibrated model scores, not measured Kd. For PEP they are not measurements of substrate binding or Km.", "",
        "| Amino acid at 177 | Variant | PEP score | Glyphosate score | Δ PEP vs WT | Δ glyphosate vs WT | Seeds (PEP/glyphosate) |",
        "|---|---|---:|---:|---:|---:|---:|"]
    def number(value):
        return "—" if value is None else f"{value:.3f}"
    def score(item):
        return "Pending" if item["n"] != 3 else f"{item['mean']:.3f} ± {item['sd']:.3f}"
    for row in rows:
        p, g = row["pep"], row["glyphosate"]
        lines.append(f"| {row['amino_acid']} — {row['name']} | {row['variant']} | {score(p)} | {score(g)} | {number(p['delta_vs_wt'])} | {number(g['delta_vs_wt'])} | {p['n']}/3; {g['n']}/3 |")
    lines += ["", "Positive Δ means weaker predicted binding than WT; negative Δ means stronger. No calibrated retention/weakening thresholds were established, so no mutation is labelled validated or selected.", ""]
    lines += ["- " + x for x in manifest["limitations"]]
    (root / "binding_table.md").write_text("\n".join(lines) + "\n")
    return report


def run(root, cache, executable):
    manifest = json.loads((root / "manifest.json").read_text())
    if importlib.metadata.version("boltz") != manifest["protocol"]["boltz"]:
        raise ValueError("Boltz version mismatch")
    if (root / "execution.json").exists():
        raise FileExistsError("Execution already started; refusing duplicate launch")
    for row in manifest["records"]:
        if sha(root / row["request"]) != row["request_sha256"] or sha(root / row["msa"]) != row["msa_sha256"]:
            raise ValueError("Registered input changed")
    events = []
    save(root / "execution.json", events)
    failures = []
    for seed in SEEDS:
        out = root / "boltz" / f"seed-{seed}"
        if out.exists():
            raise FileExistsError("Refusing previously existing predictions")
        command = [str(executable), "predict", str(root / "inputs" / f"seed-{seed}"),
            "--out_dir", str(out), "--cache", str(cache), "--model", "boltz2", "--accelerator", "gpu",
            "--devices", "1", "--seed", str(seed), "--recycling_steps", "3", "--sampling_steps", "100",
            "--diffusion_samples", "1", "--sampling_steps_affinity", "200", "--diffusion_samples_affinity", "3",
            "--num_workers", "0", "--no_kernels", "--output_format", "mmcif"]
        event = {"seed": seed, "started_at": now(), "command": command, "status": "RUNNING"}
        events.append(event)
        save(root / "execution.json", events)
        save(root / "status.json", {"state": "RUNNING", "seed": seed, "updated_at": now()})
        started = time.monotonic()
        with (root / f"seed-{seed}.log").open("w") as log:
            proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        event.update(finished_at=now(), elapsed_seconds=time.monotonic() - started,
                     returncode=proc.returncode, status="COMPLETED" if proc.returncode == 0 else "FAILED")
        save(root / "execution.json", events)
        report = collect(root)
        if proc.returncode or len([x for x in report["replicates"] if x["seed"] == seed]) != 40:
            failures.append(seed)
        print(f"Seed {seed}: {report['completed']}/120 validated results; elapsed {event['elapsed_seconds']:.1f}s", flush=True)
    report = collect(root)
    save(root / "status.json", {"state": "COMPLETED" if report["complete"] and not failures else "INCOMPLETE",
         "updated_at": now(), "completed": report["completed"], "expected": 120, "failed_seeds": failures})
    if not report["complete"] or failures:
        raise RuntimeError("Campaign incomplete; see status.json and logs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run", "collect"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--donor", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--boltz", type=Path)
    args = parser.parse_args()
    if args.mode == 'prepare' and args.donor is None:
        parser.error('prepare requires --donor')
    if args.mode == 'run' and (args.cache is None or args.boltz is None):
        parser.error('run requires --cache and --boltz')
    if args.mode == "prepare":
        prepare(args.root.resolve(), args.donor.resolve())
    elif args.mode == "run":
        run(args.root.resolve(), args.cache.resolve(), args.boltz.resolve())
    else:
        result = collect(args.root.resolve())
        print(json.dumps({k: result[k] for k in ("complete", "completed", "expected", "errors")}))
