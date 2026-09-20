"""Fresh EPSPS baseline/control calculations, with explicit model chemical states.

Structure-only predictions preserve the requested protonation hypothesis.
Affinity requests are separate because Boltz standardizes its scored ligand.
Existing MSA/reference inputs are reused and hashed; predictions are never reused.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.metadata
import itertools
import json
import math
import shutil
import subprocess
from pathlib import Path

from ..app.backends.calibration import sha, write, validate_a3m, ligand_contacts
from ..app.orchestrator.execution import ExecutionJournal, utc_now, hash_artifacts, validate_journal

SEEDS = (211, 223, 227)
SELECTED_STATES = {
    "pep": "C=C(OP(=O)([O-])[O-])C(=O)[O-]",
    "s3p": "O=C([O-])C1=C[C@@H](OP(=O)([O-])[O-])[C@@H](O)[C@H](O)C1",
    "glyphosate": "O=C([O-])C[NH2+]CP(=O)([O-])[O-]",
}
CONTEXT_LIGANDS = {
    "apo": [], "s3p_only": [("B", "s3p")],
    "herbicide": [("B", "glyphosate"), ("C", "s3p")],
    "native": [("B", "pep"), ("C", "s3p")],
    "native_s3p": [("B", "pep"), ("C", "s3p")],
}


def validate_campaign_inputs(root, manifest):
    """Validate current bytes against prospectively saved hashes; no protocol edits."""
    root = Path(root)
    checked = []
    for name, subject in manifest["subjects"].items():
        path = root / subject["msa"]["path"]
        if sha(path) != subject["msa"]["sha256"]:
            raise ValueError(f"Prospective MSA bytes changed for {name}")
        validate_a3m(path.read_text(), subject["sequence"])
        checked.append(path)
    identities = set()
    for row in manifest["records"]:
        if row["id"] in identities:
            raise ValueError("Duplicate prospective model record identity")
        identities.add(row["id"])
        request = root / row["request"]
        if sha(request) != row["request_sha256"]:
            raise ValueError("Prospective model input changed")
        data = json.loads(request.read_text())
        protein = data["sequences"][0]["protein"]
        subject = manifest["subjects"][row["subject"]]
        if protein["sequence"] != subject["sequence"] or not Path(protein["msa"]).as_posix().endswith(Path(subject["msa"]["path"]).as_posix()):
            raise ValueError("Request sequence or MSA source differs from the registered subject")
        checked.append(request)
    return {"status": "INPUT_BYTES_MATCH_REGISTERED_HASHES", "checked_artifacts": hash_artifacts(sorted(set(checked))),
            "record_count": len(identities), "scope": "Current byte integrity; no assertion about past modification times"}


def prepare(root, donor, reference, *, remote_root=None, chemical_audit=None):
    from rdkit import Chem
    root, donor = Path(root), Path(donor)
    if root.exists():
        raise FileExistsError("Fresh campaign requires a new directory")
    root.mkdir(parents=True)
    journal = ExecutionJournal(root / "preparation_journal.json")
    source = json.loads((donor / "calibration_inputs.json").read_text())
    manifest = {"schema_version": "fresh-epsps-1", "registered_at": utc_now(),
                "target_agi": "AT2G45300", "sequence": source["sequence"],
                "subjects": source["subjects"], "ligands": SELECTED_STATES, "seeds": list(SEEDS),
                "records": [], "source_manifest_sha256": sha(donor / "calibration_inputs.json"),
                "prediction_reuse": False, "ph": 7.4,
                "chemical_state_status": "declared protonation hypothesis, not measured bound-state populations",
                "protocol": {"model": "boltz2", "version": "2.2.1", "recycling_steps": 3,
                             "sampling_steps": 100, "diffusion_samples": 1,
                             "sampling_steps_affinity": 200, "diffusion_samples_affinity": 3},
                "limitations": ["No templates supplied; training-set overlap is unknown",
                    "Affinity standardization can change the scored ligand state; separate structure-only models preserve requested states",
                    "Replicate ranges are sampling variability, not calibrated confidence intervals",
                    "Apo7PXY and closed complexes represent different biological conformations",
                    "Explicit solvent, assay buffer ions and pH-dependent protein side-chain ensembles are not simulated"],
                "prospective_rules": {"heavy_atom_contact_cutoff_angstrom": 5.0,
                    "contact_reproducibility_min": 2/3, "minimum_independent_methods": 2,
                    "ligand_pose_cluster_rmsd_angstrom": 2.0,
                    "nomination_requires_validated_control_panel": True,
                    "no_candidate_derived_thresholds": True}}

    def inputs():
        (root / "msa").mkdir()
        for name, subject in manifest["subjects"].items():
            path = donor / subject["msa"]["path"]
            if sha(path) != subject["msa"]["sha256"]:
                raise ValueError("Source MSA hash mismatch")
            validate_a3m(path.read_text(), subject["sequence"])
            shutil.copy2(path, root / subject["msa"]["path"])
        shutil.copytree(donor / "references", root / "references")
        shutil.copy2(reference, root / "references/7PXY.cif")
        return {"subjects": list(manifest["subjects"]), "msa_reused_as_input_only": True}

    journal.execute("registry", inputs, inputs={"target": "AT2G45300"},
                    input_artifacts=[donor / "calibration_inputs.json", reference],
                    output_artifacts=lambda _: list((root / "msa").glob("*")) + list((root / "references").glob("*")))
    molecules = {}
    for name, smiles in SELECTED_STATES.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid ligand state")
        molecules[name] = {"smiles": Chem.MolToSmiles(mol), "formal_charge": Chem.GetFormalCharge(mol),
                           "stereo": Chem.FindMolChiralCenters(mol, includeUnassigned=True)}
    write(root / "selected_chemical_states.json", {"ligands": molecules, "ph_assumption": 7.4,
         "selection": "carboxylate and dianionic phosphate/phosphonate; protonated glyphosate amine",
         "status": manifest["chemical_state_status"]})
    if chemical_audit:
        shutil.copy2(chemical_audit, root / "chemical_state_enumeration.json")
    base = Path(remote_root) if remote_root else root.resolve()
    for mode, contexts in (("structure", ("apo", "s3p_only", "herbicide", "native")),
                           ("affinity", ("herbicide", "native", "native_s3p"))):
        for seed, (subject_name, subject), context in itertools.product(SEEDS, manifest["subjects"].items(), contexts):
            record_id = f"{subject_name}-{context}-{mode}-{seed}"
            sequences = [{"protein": {"id": "A", "sequence": subject["sequence"],
                           "msa": str(base / subject["msa"]["path"])}}]
            ligands = [{"id": chain, "name": name, "smiles": SELECTED_STATES[name]}
                       for chain, name in CONTEXT_LIGANDS[context]]
            sequences += [{"ligand": {"id": x["id"], "smiles": x["smiles"]}} for x in ligands]
            request = {"version": 1, "sequences": sequences}
            binder = "C" if context == "native_s3p" else "B"
            if mode == "affinity":
                request["properties"] = [{"affinity": {"binder": binder}}]
            path = root / "inputs" / mode / f"seed-{seed}" / f"{record_id}.yaml"
            write(path, request)
            manifest["records"].append({"id": record_id, "subject": subject_name, "context": context,
                "mode": mode, "seed": seed, "ligand": next((x["name"] for x in ligands if x["id"] == binder), None),
                "ligands": ligands, "request": str(path.relative_to(root)), "request_sha256": sha(path)})
    write(root / "calibration_inputs.json", manifest)
    journal.finalize()
    return manifest


def _run_predictions(root, cache, boltz, journal):
    from boltz.data.parse.schema import standardize
    from rdkit import Chem
    root = Path(root).resolve()
    manifest = json.loads((root / "calibration_inputs.json").read_text())
    if importlib.metadata.version("boltz") != manifest["protocol"]["version"]:
        raise ValueError("Boltz version changed")
    journal.execute("model_input_integrity", lambda: validate_campaign_inputs(root, manifest),
                    inputs={"registered_at": manifest["registered_at"]},
                    input_artifacts=[root / "calibration_inputs.json"], evidence_status="VALIDATED")
    # Record the parser's effective chemistry before any predictions execute.
    for row in manifest["records"]:
        request = root / row["request"]
        if sha(request) != row["request_sha256"]:
            raise ValueError("Prospective model input changed")
        for ligand in row["ligands"]:
            scored = row["mode"] == "affinity" and ligand["name"] == row["ligand"]
            actual = standardize(ligand["smiles"]) if scored else ligand["smiles"]
            ligand["model_smiles"] = actual
            ligand["formal_charge"] = Chem.GetFormalCharge(Chem.MolFromSmiles(actual))
            ligand["standardized_by_affinity_parser"] = scored
            ligand["requested_state_preserved"] = Chem.MolToSmiles(Chem.MolFromSmiles(actual)) == Chem.MolToSmiles(Chem.MolFromSmiles(ligand["smiles"]))
    write(root / "effective_model_inputs.json", manifest)
    failures = []
    for mode, seed in itertools.product(("structure", "affinity"), SEEDS):
        rows = [r for r in manifest["records"] if r["mode"] == mode and r["seed"] == seed]
        out = root / "boltz" / mode / f"seed-{seed}"
        if out.exists():
            raise FileExistsError("Refusing to reuse model predictions")
        command = [str(boltz), "predict", str(root / "inputs" / mode / f"seed-{seed}"),
            "--out_dir", str(out), "--cache", str(cache), "--model", "boltz2", "--accelerator", "gpu",
            "--devices", "1", "--seed", str(seed), "--recycling_steps", "3", "--sampling_steps", "100",
            "--diffusion_samples", "1", "--sampling_steps_affinity", "200", "--diffusion_samples_affinity", "3",
            "--num_workers", "0", "--no_kernels"]
        log = root / f"{mode}-{seed}.log"
        def calculate():
            with log.open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, timeout=14400)
            if result.returncode:
                raise RuntimeError(f"Boltz {mode} seed {seed} failed; see log")
            for row in rows:
                prefix = out / f"boltz_results_seed-{seed}" / "predictions" / row["id"]
                row["structure"] = str((prefix / f"{row['id']}_model_0.cif").relative_to(root))
                if not (root / row["structure"]).is_file():
                    raise FileNotFoundError("Boltz returned success but expected model is missing")
            return {"records": [r["id"] for r in rows], "exit_code": result.returncode}
        try:
            journal.execute(mode, calculate, inputs={"command": command},
                input_artifacts=[root / r["request"] for r in rows],
                output_artifacts=lambda _: [log] + [path for row in rows for path in sorted((root / row["structure"]).parent.glob("*")) if path.is_file()],
                tools=[{"name": "boltz", "version": importlib.metadata.version("boltz"), "model": "boltz2", "executable": str(boltz)}])
        except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as exc:
            failures.append({"mode": mode, "seed": seed, "error_type": type(exc).__name__, "log": str(log)})
        write(root / "effective_model_inputs.json", manifest)
    if failures:
        write(root / "model_execution_failures.json", failures)
        raise RuntimeError("One or more model batches failed; independent batches were still attempted")
    return manifest


def run(root, cache, boltz):
    root = Path(root).resolve()
    journal = ExecutionJournal(root / "model_execution_journal.json")
    try:
        return journal.execute("model_campaign", lambda: _run_predictions(root, cache, boltz, journal),
                               inputs={"root": str(root), "cache": str(cache), "boltz_executable": str(boltz)},
                               input_artifacts=[root / "calibration_inputs.json"])
    finally:
        journal.finalize()


def validate_ligand_composition(structure, record):
    """Validate component completeness against effective, not requested, chemistry.

    Atom inventory cannot prove protonation, bond order or stereochemistry;
    those remain explicit protocol hypotheses even after this check passes.
    """
    from rdkit import Chem
    expected_chains = {ligand["id"] for ligand in record["ligands"]}
    if len(expected_chains) != len(record["ligands"]) or "A" in expected_chains:
        raise ValueError("Ligand chains must be distinct from the protein chain")
    observed_chains = {chain.id for chain in structure[0]} - {"A"}
    if observed_chains != expected_chains:
        raise ValueError("Predicted ligand chains do not match the registered complex (including apo absence)")
    inventories = {}
    for ligand in record["ligands"]:
        if not ligand.get("model_smiles"):
            raise ValueError("Effective model SMILES required for ligand atom validation")
        molecule = Chem.MolFromSmiles(ligand["model_smiles"])
        if molecule is None:
            raise ValueError("Invalid effective model ligand state")
        expected = Counter(atom.GetSymbol().upper() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1)
        atoms = list(structure[0][ligand["id"]].get_atoms())
        observed = Counter()
        for atom in atoms:
            element = str(atom.element or "").upper().strip()
            if not element or element == "X":
                raise ValueError("Predicted ligand atom has unknown chemical element")
            if element not in ("H", "D"):
                observed[element] += 1
        if sum(observed.values()) != molecule.GetNumHeavyAtoms() or observed != expected:
            raise ValueError(f"Ligand heavy-atom inventory does not match effective SMILES for {record['id']}/{ligand['id']}")
        inventories[ligand["id"]] = {"ligand": ligand["name"], "expected_heavy_atoms": molecule.GetNumHeavyAtoms(),
                                     "observed_heavy_atoms": sum(observed.values()), "expected_elements": dict(sorted(expected.items())),
                                     "observed_elements": dict(sorted(observed.items())), "effective_model_smiles": ligand["model_smiles"],
                                     "validation_scope": "heavy-atom and element completeness only; no stereochemical/protonation inference"}
    return inventories


def verify_originating_model_hashes(root, manifest):
    """Match every downloaded CIF to exactly one actual originating model event."""
    root = Path(root).resolve()
    path = root / "model_execution_journal.json"
    original = json.loads(path.read_text())
    validate_journal(original)
    source_artifacts = [(event, artifact) for event in original["stages"]
                        if event["execution_status"] == "COMPLETED" and event["execution_mode"] == "executed"
                        for artifact in event.get("output_artifacts", []) if artifact["path"].endswith(".cif")]
    if len(source_artifacts) != len(manifest["records"]):
        raise ValueError("Originating journal must contain exactly the registered cohort of CIF outputs")
    validated, matched, seen = [], set(), set()
    for row in manifest["records"]:
        relative = Path(row["structure"])
        if relative.is_absolute() or ".." in relative.parts or row["structure"] in seen:
            raise ValueError("Expected unique campaign-relative CIF paths")
        seen.add(row["structure"])
        suffix = "/" + relative.as_posix()
        matches = [(event, artifact) for event, artifact in source_artifacts if artifact["path"].endswith(suffix)]
        if len(matches) != 1:
            raise ValueError("Each downloaded CIF needs exactly one unambiguous originating event hash")
        event, expected = matches[0]
        actual = hash_artifacts([root / relative])[0]
        if actual["sha256"] != expected["sha256"] or actual["size_bytes"] != expected["size_bytes"]:
            raise ValueError(f"Downloaded CIF differs from originating model bytes: {row['id']}")
        identity = (event["id"], expected["path"])
        if identity in matched:
            raise ValueError("One originating model artifact cannot validate two output records")
        matched.add(identity)
        validated.append({"record_id": row["id"], "structure": relative.as_posix(), "source_event_id": event["id"],
                          "source_run_id": original["run_id"], "sha256": actual["sha256"], "size_bytes": actual["size_bytes"]})
    return {"status": "ALL_CIFS_MATCH_ORIGINATING_EXECUTION_HASHES", "verified_cif_count": len(validated),
            "originating_journal_sha256": sha(path), "records": validated, "originating_timestamps_modified": False}


def _analyze_predictions(root, manifest):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    from .epsps_local_campaign import validate_affinity
    root = Path(root)
    predictions = []
    for row in manifest["records"]:
        structure_path = root / row["structure"]
        structure = MMCIFParser(QUIET=True).get_structure(row["id"], structure_path)
        expected = manifest["subjects"][row["subject"]]["sequence"]
        actual = "".join(seq1(r.resname) for r in structure[0]["A"] if r.id[0] == " ")
        if actual != expected or any(not math.isfinite(float(x)) for a in structure.get_atoms() for x in a.coord):
            raise ValueError("Output structure failed sequence/finite-coordinate validation")
        ligand_inventory = validate_ligand_composition(structure, row)
        contacts = {x["id"]: ligand_contacts(structure, "A", [x["id"]]) for x in row["ligands"]}
        prediction = {**row, "structure_sha256": sha(structure_path), "contacts": contacts, "ligand_atom_validation": ligand_inventory}
        if row["mode"] == "affinity":
            path = structure_path.parent / f"affinity_{row['id']}.json"
            prediction["predicted_pIC50"] = validate_affinity(json.loads(path.read_text()))
            prediction["affinity_sha256"] = sha(path)
        predictions.append(prediction)
    write(root / "predictions.json", predictions)
    return predictions


def _analyze_postflight(root, journal):
    """New postflight validation event; never reconstruct missing model runtimes."""
    root = Path(root).resolve()
    manifest_path = root / "effective_model_inputs.json"
    manifest = json.loads(manifest_path.read_text())
    registered = json.loads((root / "calibration_inputs.json").read_text())
    # Effective-input annotations are additive. Identity/protocol cannot change.
    for key in ("registered_at", "subjects", "seeds", "protocol", "sequence"):
        if manifest.get(key) != registered.get(key):
            raise ValueError(f"Effective output manifest changed registered {key}")
    journal.execute("postflight_input_integrity", lambda: validate_campaign_inputs(root, registered),
                    inputs={"registered_at": registered["registered_at"], "scope": "postflight byte verification"},
                    input_artifacts=[root / "calibration_inputs.json", manifest_path], evidence_status="VALIDATED")
    registered_rows = {row["id"]: row for row in registered["records"]}
    observed = [row["id"] for row in manifest["records"]]
    if len(observed) != len(set(observed)) or set(observed) != set(registered_rows):
        raise ValueError("Output manifest does not cover the complete registered model cohort")
    journal.execute("postflight_origin_hash_verification", lambda: verify_originating_model_hashes(root, manifest),
                    inputs={"registered_records": observed}, input_artifacts=[root / "model_execution_journal.json"] + [root / row["structure"] for row in manifest["records"]],
                    evidence_status="VALIDATED")
    source_files = [manifest_path, root / "calibration_inputs.json", root / "model_execution_journal.json"]
    source_files.extend(root / subject["msa"]["path"] for subject in manifest["subjects"].values())
    for row in manifest["records"]:
        planned = registered_rows[row["id"]]
        for key in ("subject", "context", "mode", "seed", "request", "request_sha256", "ligand"):
            if row.get(key) != planned.get(key):
                raise ValueError(f"Output record changed registered {key}")
        if len(row["ligands"]) != len(planned["ligands"]) or any(
                any(actual.get(key) != original.get(key) for key in ("id", "name", "smiles"))
                for actual, original in zip(row["ligands"], planned["ligands"])):
            raise ValueError("Effective-input annotations changed a registered ligand identity/state")
        structure = root / row["structure"]
        source_files.extend([root / row["request"], structure])
        confidence = structure.parent / f"confidence_{row['id']}_model_0.json"
        if not confidence.is_file():
            raise FileNotFoundError("Expected model-confidence output missing")
        source_files.append(confidence)
        if row["mode"] == "affinity":
            source_files.append(structure.parent / f"affinity_{row['id']}.json")
    tools = []
    for name in ("biopython", "numpy"):
        try:
            tools.append({"name": name, "version": importlib.metadata.version(name), "scope": "postflight analysis environment"})
        except importlib.metadata.PackageNotFoundError:
            pass
    return journal.execute("model_output_validation", lambda: _analyze_predictions(root, manifest),
                           inputs={"scope": "freshly executed validation of existing model outputs; not model re-execution", "records": observed},
                           input_artifacts=sorted(set(source_files)), output_artifacts=[root / "predictions.json"], tools=tools,
                           evidence_status="AVAILABLE")


def analyze(root):
    """Record all postflight failures without changing originating model events."""
    root = Path(root).resolve()
    journal = ExecutionJournal(root / "postflight_execution_journal.json")
    try:
        return journal.execute("postflight_analysis", lambda: _analyze_postflight(root, journal),
                               inputs={"scope": "current artifact integrity and output validation"},
                               input_artifacts=[root / "calibration_inputs.json", root / "effective_model_inputs.json"],
                               output_artifacts=[root / "predictions.json"])
    finally:
        journal.finalize()


def verify_downloaded_outputs(root):
    """Add independent origin/composition verification without rewriting old evidence."""
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    from rdkit import rdBase
    root = Path(root).resolve()
    manifest = json.loads((root / "effective_model_inputs.json").read_text())
    report_path = root / "source_integrity_report.json"
    journal = ExecutionJournal(root / "source_integrity_verification.json")
    originals = {name: sha(root / name) for name in ("model_execution_journal.json", "postflight_execution_journal.json", "predictions.json")}
    def verify():
        origin = verify_originating_model_hashes(root, manifest)
        parser = MMCIFParser(QUIET=True)
        rows = []
        for row in manifest["records"]:
            structure = parser.get_structure(row["id"], root / row["structure"])
            actual = "".join(seq1(r.resname) for r in structure[0]["A"] if r.id[0] == " ")
            if actual != manifest["subjects"][row["subject"]]["sequence"] or any(not math.isfinite(float(x)) for a in structure.get_atoms() for x in a.coord):
                raise ValueError("Downloaded protein failed sequence or finite-coordinate validation")
            rows.append({"record_id": row["id"], "context": row["context"], "mode": row["mode"],
                         "ligand_inventory": validate_ligand_composition(structure, row)})
        if any(sha(root / name) != expected for name, expected in originals.items()):
            raise RuntimeError("Original model/postflight artifacts changed during verification")
        result = {"status": "ORIGINATING_CIF_HASHES_AND_LIGAND_INVENTORIES_VERIFIED", "cohort_count": len(rows),
                  "origin_hash_validation": origin, "ligand_atom_validation": rows,
                  "preserved_original_artifact_hashes": originals,
                  "scope": "Fresh integrity/atom-count verification only; no model rerun or edited originating timestamps",
                  "limitations": ["Atom counts and elemental composition do not establish stereochemistry, bond order, protonation state, affinity or enzyme function",
                                  "Original model journal contains CIF hashes; other output artifacts retain separately measured postflight hashes"]}
        write(report_path, result)
        return result
    try:
        return journal.execute("model_output_integrity", verify,
            inputs={"purpose": "verify original model output bytes and effective ligand atom completeness", "expected_records": len(manifest["records"])},
            input_artifacts=[root / "effective_model_inputs.json", root / "model_execution_journal.json", root / "postflight_execution_journal.json", root / "predictions.json"] + [root / row["structure"] for row in manifest["records"]],
            output_artifacts=[report_path], evidence_status="VALIDATED",
            tools=[{"name": "biopython", "version": importlib.metadata.version("biopython")}, {"name": "rdkit", "version": rdBase.rdkitVersion}])
    finally:
        journal.finalize()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "run", "analyze", "verify"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--donor", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--remote-root")
    parser.add_argument("--chemical-audit", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--boltz", default="boltz")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.root, args.donor, args.reference, remote_root=args.remote_root, chemical_audit=args.chemical_audit)
    elif args.action == "run":
        run(args.root, args.cache, args.boltz)
    elif args.action == "analyze":
        analyze(args.root)
    else:
        verify_downloaded_outputs(args.root)


if __name__ == "__main__":
    main()
