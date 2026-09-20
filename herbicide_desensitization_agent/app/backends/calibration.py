"""Reproducible EPSPS WT/experimental-control calibration, separate from design."""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import statistics
import subprocess
import time
import shutil
from pathlib import Path
from urllib.request import urlopen

from .structural import TMAlignStructuralMatcher


LEGEND = (
    "Calibration uses Arabidopsis EPSPS WT and experimental E. coli homolog controls, not new mutation designs. "
    "Three seeded predictions per condition estimate sampling variability, not calibrated confidence intervals. "
    "TM-scores and CA lDDT range from 0 to 1 (higher is closer); globally fitted CA and pocket-backbone RMSD "
    "are in angstroms (lower is closer). All replicate pairs are retained. Experimental comparisons use the "
    "same observed protein sequence and ligand context as the reference crystal. E. coli control accuracy "
    "does not establish Arabidopsis accuracy. Predicted pIC50 is 6 minus Boltz affinity_pred_value; higher "
    "means stronger predicted binding. It is not measured Kd, Km, Ki or enzyme activity. Templates are not "
    "supplied or forced. DiffDock is diagnostic only: its protein-only receptor cannot establish agreement "
    "for S3P-containing complexes. Structural and functional control outcomes are reported separately."
    " Training-set exclusion has not been established for these historical experimental controls."
)
SEEDS = (101, 103, 107)
CONTEXTS = {"arabidopsis_wt": ("herbicide", "native", "native_s3p"),
            "ecoli_wt": ("herbicide", "native"), "ecoli_G96A": ("herbicide", "native", "s3p_only")}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def validate_record_coverage(manifest):
    if tuple(manifest["seeds"]) != SEEDS or set(manifest["subjects"]) != set(CONTEXTS):
        raise ValueError("Calibration subject/seed protocol changed")
    expected = {(name, context, seed) for name, contexts in CONTEXTS.items() for context in contexts for seed in SEEDS}
    rows = manifest["records"]
    actual = {(r["subject"], r["context"], r["seed"]) for r in rows}
    if actual != expected or len(rows) != len(expected) or any(r["id"] != f"{r['subject']}-{r['context']}-{r['seed']}" for r in rows):
        raise ValueError("Missing, duplicate or mislabeled calibration records")


def validate_a3m(text, sequence):
    from io import StringIO
    from Bio import SeqIO

    records = list(SeqIO.parse(StringIO(text), "fasta"))
    aligned = ["".join(c for c in str(r.seq) if not c.islower() and c != ".") for r in records]
    if not aligned or aligned[0] != sequence:
        raise ValueError("MSA query does not exactly match modeled sequence")
    if any(len(row) != len(sequence) or set(row) - set("ACDEFGHIKLMNPQRSTVWYX-") for row in aligned):
        raise ValueError("Invalid A3M alignment columns or residues")
    if len(set(aligned)) < 2:
        raise ValueError("MSA has no distinct homologs; query-only fallback is prohibited")
    return {"rows": len(records), "unique_aligned_rows": len(set(aligned)), "query_length": len(sequence)}


def replace_msa_query(text, old_sequence, new_sequence):
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq

    validate_a3m(text, old_sequence)
    if len(old_sequence) != len(new_sequence):
        raise ValueError("Matched control MSA requires equal sequence lengths")
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    records[0].seq = Seq(new_sequence)
    output = StringIO()
    SeqIO.write(records, output, "fasta-2line")
    result = output.getvalue()
    validate_a3m(result, new_sequence)
    return result


def ligand_contacts(structure, chain, ligand_chains, cutoff=5.0):
    import numpy as np

    ligand_atoms = [a.coord for c in ligand_chains for a in structure[0][c].get_atoms() if a.element not in {"H", "D"}]
    if not ligand_atoms:
        raise ValueError("Missing ligand heavy atoms")
    xyz = np.asarray(ligand_atoms)
    return [r.id[1] for r in structure[0][chain] if r.id[0] == " " and
            any(np.linalg.norm(xyz - a.coord, axis=1).min() <= cutoff for a in r if a.element not in {"H", "D"})]


class EPSPSCalibrationRunner:
    def __init__(self, previous, output, boltz_cache, boltz_executable="boltz", reuse_calibration=None):
        self.previous, self.root, self.cache = Path(previous), Path(output), Path(boltz_cache)
        self.boltz_executable = str(boltz_executable)
        self.reuse_calibration = Path(reuse_calibration) if reuse_calibration else None

    def prepare(self):
        from Bio.PDB import MMCIFParser, PDBIO
        from Bio.SeqUtils import seq1
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        if self.root.exists():
            raise ValueError("Calibration requires a fresh output directory")
        previous = json.loads((self.previous / "inputs/input_manifest.json").read_text())
        if (previous.get("target_agi") != "AT2G45300" or len(previous.get("sequence", "")) != 520 or
                set(previous["sequence"]) - set("ACDEFGHIKLMNPQRSTVWY") or
                set(previous.get("ligands", {})) != {"pep", "s3p", "glyphosate"}):
            raise ValueError("Expected the validated EPSPS sequence and all three ligands")
        self.root.mkdir(parents=True)
        controls = {}
        for code, label in (("1G6S", "ecoli_wt"), ("1MI4", "ecoli_G96A")):
            url = f"https://files.rcsb.org/download/{code}.cif"
            source = self.root / "references" / f"{code}.cif"
            source.parent.mkdir(exist_ok=True)
            with urlopen(url, timeout=90) as response:
                source.write_bytes(response.read())
            model = MMCIFParser(QUIET=True).get_structure(code, source)
            residues = [r for r in model[0]["A"] if r.id[0] == " " and "CA" in r]
            sequence = "".join(seq1(r.resname) for r in residues)
            original_numbers = [r.id[1] for r in residues]
            # Normalize observed-residue numbering explicitly for sequence-matched comparisons.
            protein = model[0]["A"].copy()
            for r in list(protein):
                protein.detach_child(r.id)
            for index, r in enumerate(residues, 1):
                copy = r.copy()
                copy.id = (" ", index, " ")
                protein.add(copy)
            normalized = self.root / "references" / f"{code}-observed-protein.pdb"
            writer = PDBIO()
            writer.set_structure(protein)
            writer.save(str(normalized))
            ligand_xyz = [a.coord for r in model[0]["A"] if r.resname in {"S3P", "GPJ"} for a in r
                          if a.element not in {"H", "D"}]
            import numpy as np
            if not ligand_xyz:
                raise ValueError("Reference ligand context is absent")
            pocket = [i for i, r in enumerate(residues, 1)
                      if any(np.linalg.norm(np.asarray(ligand_xyz) - a.coord, axis=1).min() <= 5.0 for a in r
                             if a.element not in {"H", "D"})]
            controls[label] = {"sequence": sequence, "reference": str(normalized.relative_to(self.root)),
                               "pdb": code, "source_url": url, "reference_sha256": sha(normalized),
                               "source_sha256": sha(source), "observed_author_residue_numbers": original_numbers,
                               "pocket_residues": pocket, "scope": "E. coli homolog method control; not Arabidopsis evidence"}
        wt, mutant = controls["ecoli_wt"]["sequence"], controls["ecoli_G96A"]["sequence"]
        changes = [(i, a, b) for i, (a, b) in enumerate(zip(wt, mutant), 1) if a != b]
        if len(wt) != len(mutant) or len(changes) != 1 or changes[0][1:] != ("G", "A"):
            raise ValueError("Experimental control sequences are not a matched G-to-A pair")
        full = previous["sequence"]
        subjects = {"arabidopsis_wt": {"sequence": full[76:], "scope": "Arabidopsis WT", "offset": 76}, **controls}
        chemical_audit = {}
        for name, smiles in previous["ligands"].items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid chemical input")
            chemical_audit[name] = {"smiles": Chem.MolToSmiles(mol), "formal_charge": Chem.GetFormalCharge(mol),
                                    "formula": rdMolDescriptors.CalcMolFormula(mol), "molecular_weight": Descriptors.MolWt(mol),
                                    "stereocentres": Chem.FindMolChiralCenters(mol, includeUnassigned=True),
                                    "chemical_identity_check": "RDKit parse and stereo audit",
                                    "ph": None, "protonation_validated": False,
                                    "limitation": "Stored PubChem chemical state; no justified assay pH or validated protonation ensemble supplied"}
        manifest = {"target_agi": "AT2G45300", "sequence": full, "subjects": subjects,
                    "ligands": previous["ligands"], "seeds": list(SEEDS), "controls": controls,
                    "chemical_audit": chemical_audit, "legend": LEGEND,
                    "diffdock_context": "protein-only diagnostic; not a matched S3P-aware oracle"}
        manifest["source_input_manifest_sha256"] = sha(self.previous / "inputs/input_manifest.json")
        write(self.root / "calibration_inputs.json", manifest)
        write(self.root / "chemical_state_audit.json", {"ligands": chemical_audit, "legend": LEGEND})
        sources = [{"id": c["pdb"], "url": "https://www.rcsb.org/structure/" + c["pdb"],
                    "artifact": f"references/{c['pdb']}.cif", "sha256": c["source_sha256"], "scope": c["scope"]}
                   for c in controls.values()]
        write(self.root / "evidence.json", {"target_agi": "AT2G45300", "sources": sources, "claims": [
            {"text": "1G6S is an experimental E. coli EPSPS structure containing S3P and glyphosate; it is a homolog control, not an Arabidopsis reference.", "source_ids": ["1G6S"]},
            {"text": "1MI4 is an experimental E. coli G96A EPSPS structure containing S3P. Its associated study reports glyphosate insensitivity with impaired PEP affinity; kinetics are not Kd labels.", "source_ids": ["1MI4"]}],
            "limitations": ["This historical protocol does not include the available same-target open-apo structure 7PXY or a complete Arabidopsis functional control panel. Use the fresh campaign for state-matched apo comparisons."]})
        return manifest

    def prepare_msa(self):
        manifest = json.loads((self.root / "calibration_inputs.json").read_text())
        subjects = manifest["subjects"]
        names = ["arabidopsis_wt", "ecoli_wt"]
        if self.reuse_calibration:
            donor = json.loads((self.reuse_calibration / "calibration_inputs.json").read_text())
            donor_report = json.loads((self.reuse_calibration / "calibration_report.json").read_text())
            for artifact in donor_report["artifacts"]:
                if sha(self.reuse_calibration / artifact["path"]) != artifact["sha256"]:
                    raise ValueError("Reuse source artifact changed")
            texts = []
            for name in names:
                if donor["subjects"][name]["sequence"] != subjects[name]["sequence"]:
                    raise ValueError("Reuse source sequence mismatch")
                texts.append((self.reuse_calibration / donor["subjects"][name]["msa"]["path"]).read_text())
        else:
            from boltz.main import run_mmseqs2
            texts = run_mmseqs2([subjects[n]["sequence"] for n in names], self.root / "msa-search",
                                use_env=True, use_pairing=False, host_url="https://api.colabfold.com")
        names.append("ecoli_G96A")
        texts = list(texts) + [replace_msa_query(texts[1], subjects["ecoli_wt"]["sequence"], subjects["ecoli_G96A"]["sequence"])]
        msa_dir = self.root / "msa"
        msa_dir.mkdir(exist_ok=True)
        for name, alignment in zip(names, texts, strict=True):
            stats = validate_a3m(alignment, subjects[name]["sequence"])
            path = msa_dir / f"{name}.a3m"
            path.write_text(alignment)
            subjects[name]["msa"] = {"path": str(path.relative_to(self.root)), "sha256": sha(path), **stats,
                                      "method": "ColabFold MMseqs2 server, unpaired, environmental sequences enabled"}
            if name == "ecoli_G96A":
                subjects[name]["msa"]["method"] = "WT homolog alignment, only the query changed to the experimental mutant sequence"
                subjects[name]["msa"]["homolog_rows_shared_with"] = "ecoli_wt"
        records = []
        for name, subject in subjects.items():
            contexts = [("herbicide", "glyphosate", True), ("native", "pep", True)]
            if name == "arabidopsis_wt":
                contexts.append(("native_s3p", "pep", True))
            if name == "ecoli_G96A":
                contexts.append(("s3p_only", "s3p", False))
            for seed in manifest["seeds"]:
                for context, primary, cosubstrate in contexts:
                    record_id = f"{name}-{context}-{seed}"
                    sequences = [{"protein": {"id": "A", "sequence": subject["sequence"],
                                                "msa": str((self.root / subject["msa"]["path"]).resolve())}},
                                 {"ligand": {"id": "B", "smiles": manifest["ligands"][primary]}}]
                    if cosubstrate:
                        sequences.append({"ligand": {"id": "C", "smiles": manifest["ligands"]["s3p"]}})
                    request = {"version": 1, "sequences": sequences,
                               "properties": [{"affinity": {"binder": "C" if context == "native_s3p" else "B"}}]}
                    path = self.root / "inputs" / f"seed-{seed}" / f"{record_id}.yaml"
                    write(path, request)
                    records.append({"id": record_id, "subject": name, "context": context, "seed": seed,
                                    "ligand": "s3p" if context in {"native_s3p", "s3p_only"} else primary,
                                    "request": str(path.relative_to(self.root)), "request_sha256": sha(path)})
        manifest["records"] = records
        if self.reuse_calibration:
            donor_rows = {r["id"]: r for r in donor["records"]}
            execution = json.loads((self.reuse_calibration / "execution.json").read_text())
            expected_flags = {"--model": "boltz2", "--recycling_steps": "3", "--sampling_steps": "100",
                              "--diffusion_samples": "1", "--sampling_steps_affinity": "200", "--diffusion_samples_affinity": "3"}
            for row in records:
                if row["subject"] == "ecoli_G96A":
                    continue
                source_row = donor_rows[row["id"]]
                new_request = json.loads((self.root / row["request"]).read_text())
                old_request = json.loads((self.reuse_calibration / source_row["request"]).read_text())
                for request, base in ((new_request, self.root), (old_request, self.reuse_calibration)):
                    msa = request["sequences"][0]["protein"].pop("msa")
                    # MSA bytes, rather than machine-specific paths, define this part of the protocol.
                    request["sequences"][0]["protein"]["msa_sha256"] = sha(base / "msa" / Path(msa).name)
                if new_request != old_request or row["seed"] != source_row["seed"]:
                    raise ValueError("Cannot reuse predictions with a different request or seed")
                stage = next(s for s in execution["stages"] if s["stage"] == f"boltz-{row['seed']}")
                command = stage["command"]
                if (stage["exit_code"] != 0 or command[command.index("--seed") + 1] != str(row["seed"]) or
                        any(command[command.index(k) + 1] != v for k, v in expected_flags.items())):
                    raise ValueError("Cannot reuse predictions from another inference protocol")
                relative = Path("boltz") / f"seed-{row['seed']}" / f"boltz_results_seed-{row['seed']}" / "predictions" / row["id"]
                shutil.copytree(self.reuse_calibration / relative, self.root / relative)
                row["reuse"] = {"source": str(self.reuse_calibration / relative),
                               "source_report_sha256": sha(self.reuse_calibration / "calibration_report.json")}
            for row in records:
                if "reuse" not in row:
                    destination = self.root / "run-inputs" / f"seed-{row['seed']}" / Path(row["request"]).name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(self.root / row["request"], destination)
            manifest["reused_prediction_count"] = sum("reuse" in row for row in records)
        manifest["new_prediction_count"] = sum("reuse" not in row for row in records)
        validate_record_coverage(manifest)
        write(self.root / "calibration_inputs.json", manifest)
        return manifest

    def run(self):
        import importlib.metadata
        import importlib.util

        if importlib.metadata.version("boltz") != "2.2.1":
            raise ValueError("Calibration protocol requires boltz==2.2.1")
        missing = [name for name in ("Bio", "numpy", "rdkit", "tmtools") if importlib.util.find_spec(name) is None]
        if missing:
            raise RuntimeError("Install calibration dependencies before inference: " + ", ".join(missing))
        if (self.root / "execution.json").exists():
            raise ValueError("Refusing to reuse a launched calibration run")
        manifest = self.prepare() if not self.root.exists() else json.loads((self.root / "calibration_inputs.json").read_text())
        state = {"status": "RUNNING", "stages": [], "legend": LEGEND}
        write(self.root / "execution.json", state)
        try:
            manifest = self.prepare_msa()
            state["predictions"] = {"new": manifest["new_prediction_count"], "reused": manifest.get("reused_prediction_count", 0)}
            state["stages"].append({"stage": "msa", "status": "COMPLETED"})
            write(self.root / "execution.json", state)
            for seed in manifest["seeds"]:
                command = [self.boltz_executable, "predict", str(self.root / ("run-inputs" if self.reuse_calibration else "inputs") / f"seed-{seed}"),
                           "--out_dir", str(self.root / "boltz" / f"seed-{seed}"), "--cache", str(self.cache),
                           "--model", "boltz2", "--accelerator", "gpu", "--devices", "1", "--seed", str(seed),
                           "--recycling_steps", "3", "--sampling_steps", "100", "--diffusion_samples", "1",
                           "--sampling_steps_affinity", "200", "--diffusion_samples_affinity", "3", "--num_workers", "0", "--no_kernels"]
                row = {"stage": f"boltz-{seed}", "status": "RUNNING", "command": command, "started": time.time()}
                state["stages"].append(row)
                write(self.root / "execution.json", state)
                with (self.root / f"boltz-{seed}.log").open("w") as log:
                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
                row.update(status="COMPLETED" if result.returncode == 0 else "FAILED", exit_code=result.returncode)
                row["finished"] = time.time()
                write(self.root / "execution.json", state)
                if result.returncode:
                    raise RuntimeError("Boltz failed; see calibration log")
            self.analyze()
            state.update(json.loads((self.root / "execution.json").read_text()))
            state["status"] = "COMPUTATION_COMPLETE"
        except BaseException as exc:
            state["status"] = "FAILED"
            state["failure_type"] = type(exc).__name__
            raise
        finally:
            write(self.root / "execution.json", state)

    def analyze(self):
        from Bio.PDB import MMCIFParser
        from Bio.SeqUtils import seq1
        from rdkit import Chem
        from ...examples.epsps_local_campaign import validate_affinity

        manifest = json.loads((self.root / "calibration_inputs.json").read_text())
        validate_record_coverage(manifest)
        artifacts, predictions, comparisons = [], [], []
        matcher = TMAlignStructuralMatcher()
        for subject in manifest["subjects"].values():
            msa = self.root / subject["msa"]["path"]
            if sha(msa) != subject["msa"]["sha256"]:
                raise ValueError("MSA changed")
            validate_a3m(msa.read_text(), subject["sequence"])
            artifacts.append({"path": str(msa.relative_to(self.root)), "sha256": sha(msa)})
        from io import StringIO
        from Bio import SeqIO
        control_alignments = [list(SeqIO.parse(StringIO((self.root / manifest["subjects"][n]["msa"]["path"]).read_text()), "fasta"))
                              for n in ("ecoli_wt", "ecoli_G96A")]
        matched_control_msa = ([str(r.seq) for r in control_alignments[0][1:]] == [str(r.seq) for r in control_alignments[1][1:]])
        for row in manifest["records"]:
            subject = manifest["subjects"][row["subject"]]
            request = self.root / row["request"]
            if sha(request) != row["request_sha256"]:
                raise ValueError("Calibration request changed")
            run_request = self.root / "run-inputs" / f"seed-{row['seed']}" / request.name
            if run_request.exists():
                if sha(run_request) != row["request_sha256"]:
                    raise ValueError("Executed calibration request changed")
                artifacts.append({"path": str(run_request.relative_to(self.root)), "sha256": sha(run_request)})
            path = self.root / "boltz" / f"seed-{row['seed']}" / f"boltz_results_seed-{row['seed']}" / "predictions" / row["id"]
            cif = path / f"{row['id']}_model_0.cif"
            structure = MMCIFParser(QUIET=True).get_structure(row["id"], cif)
            observed = "".join(seq1(r.resname) for r in structure[0]["A"] if r.id[0] == " ")
            if observed != subject["sequence"]:
                raise ValueError("Calibration output sequence mismatch")
            if any(not math.isfinite(float(x)) for a in structure.get_atoms() for x in a.coord):
                raise ValueError("Non-finite calibration coordinates")
            confidence_path = path / f"confidence_{row['id']}_model_0.json"
            confidence = json.loads(confidence_path.read_text()).get("confidence_score")
            if isinstance(confidence, bool) or not isinstance(confidence, (int, float)) or not 0 <= confidence <= 1:
                raise ValueError("Invalid calibration structure confidence")
            request_data = json.loads(request.read_text())
            protein = request_data["sequences"][0]["protein"]
            binder = "C" if row["context"] == "native_s3p" else "B"
            if (protein["sequence"] != subject["sequence"] or protein.get("msa") == "empty" or
                    request_data.get("templates") or request_data.get("constraints") or
                    request_data.get("properties") != [{"affinity": {"binder": binder}}]):
                raise ValueError("Unexpected calibration sequence, template, constraint or affinity request")
            for item in request_data["sequences"][1:]:
                ligand = item["ligand"]
                count = sum(a.element not in {"H", "D"} for a in structure[0][ligand["id"]].get_atoms())
                if count != Chem.MolFromSmiles(ligand["smiles"]).GetNumHeavyAtoms():
                    raise ValueError("Calibration ligand atom coverage mismatch")
            affinity_path = path / f"affinity_{row['id']}.json"
            affinity = json.loads(affinity_path.read_text())
            contacts = {item["ligand"]["id"]: ligand_contacts(structure, "A", [item["ligand"]["id"]])
                        for item in request_data["sequences"][1:]}
            predictions.append({**row, "structure": str(cif.relative_to(self.root)), "contacts": contacts,
                                "predicted_pIC50": validate_affinity(affinity)})
            for artifact in (cif, request, affinity_path, confidence_path):
                artifacts.append({"path": str(artifact.relative_to(self.root)), "sha256": sha(artifact)})
            reference_matches = ((row["subject"] == "ecoli_wt" and row["context"] == "herbicide") or
                                 (row["subject"] == "ecoli_G96A" and row["context"] == "s3p_only"))
            if reference_matches:
                ref = self.root / subject["reference"]
                if sha(ref) != subject["reference_sha256"]:
                    raise ValueError("Experimental reference changed")
                artifacts.append({"path": subject["reference"], "sha256": sha(ref)})
                metrics = matcher.compare(cif, ref, target_sequence=subject["sequence"],
                                          active_site_residues=subject["pocket_residues"])
                comparisons.append({"kind": "experimental_accuracy_homolog", "subject": row["subject"],
                                    "context": row["context"], "query": row["id"], "reference": subject["pdb"], **metrics})
        for name, subject in manifest["subjects"].items():
            for context in sorted({r["context"] for r in predictions if r["subject"] == name}):
                rows = [r for r in predictions if r["subject"] == name and r["context"] == context]
                if sorted(r["seed"] for r in rows) != sorted(manifest["seeds"]):
                    raise ValueError("Missing or duplicated calibration seeds")
                site = sorted({p for r in rows for ps in r["contacts"].values() for p in ps})
                for a, b in itertools.combinations(rows, 2):
                    metrics = matcher.compare(self.root / a["structure"], self.root / b["structure"],
                                              target_sequence=subject["sequence"], active_site_residues=site)
                    comparisons.append({"kind": "repeatability", "subject": name, "context": context,
                                        "query": a["id"], "reference": b["id"], **metrics})
        thresholds = {"query_tm_score": .95, "target_tm_score": .95, "alignment_lddt": .90,
                      "alignment_coverage": .95, "global_ca_rmsd_angstrom": 1.0, "active_site_rmsd_angstrom": 1.0}
        for row in comparisons:
            row["failed_metrics"] = [k for k, bound in thresholds.items() if row[k] is None or
                                     (row[k] > bound if "rmsd" in k else row[k] < bound)]
        wt_pairs = [r for r in comparisons if r["kind"] == "repeatability" and r["subject"] == "arabidopsis_wt"]
        homolog_accuracy = [r for r in comparisons if r["kind"] == "experimental_accuracy_homolog"]
        controls = []
        for ligand in ("glyphosate", "pep"):
            values = {name: [r["predicted_pIC50"] for r in predictions if r["subject"] == name and r["ligand"] == ligand]
                      for name in ("ecoli_wt", "ecoli_G96A")}
            delta = statistics.mean(values["ecoli_G96A"]) - statistics.mean(values["ecoli_wt"])
            controls.append({"ligand": ligand, "mean_mutant_minus_wt_pIC50": delta,
                             "predicted_weaker_in_mutant": delta < 0, "seed_values": values,
                             "interpretation": "Direction-only diagnostic; literature kinetic endpoints are not Boltz pIC50 or Kd"})
        wt_predictions = [r for r in predictions if r["subject"] == "arabidopsis_wt"]
        herbicide_contacts = {p + 76 for r in wt_predictions if r["context"] == "herbicide" for p in r["contacts"]["B"]}
        native_contacts = {p + 76 for r in wt_predictions for chain, positions in r["contacts"].items()
                           if chain == "C" or r["context"] != "herbicide" for p in positions}
        contact_audit = {"herbicide_contacts": sorted(herbicide_contacts), "native_contacts": sorted(native_contacts),
                         "herbicide_only_positions": sorted(herbicide_contacts - native_contacts),
                         "legend": "Union of 5 A heavy-atom contacts across WT Boltz seeds, full-sequence numbering (+76). "
                                   "Herbicide-only means outside every native-ligand contact set; this is the existing conservative "
                                   "selection rule, not proof that all shared-site substitutions are biologically impossible."}
        report = {"target_agi": manifest["target_agi"],
                  "sequence_sha256": hashlib.sha256(manifest["sequence"].encode()).hexdigest(),
                  "checks": {"msa_validated": True, "no_forced_templates": True,
                             "chemistry_reviewed": False, "wt_repeatability": bool(wt_pairs) and all(not r["failed_metrics"] for r in wt_pairs),
                             "experimental_accuracy": False, "functional_controls": False},
                  "homolog_structure_accuracy_passed": bool(homolog_accuracy) and all(not r["failed_metrics"] for r in homolog_accuracy),
                  "control_msa_homolog_rows_matched": matched_control_msa,
                  "predictions": {"total": len(predictions), "new": manifest.get("new_prediction_count", len(predictions)),
                                  "reused": manifest.get("reused_prediction_count", 0)},
                  "thresholds": thresholds, "comparisons": comparisons, "control_directions": controls,
                  "contact_audit": contact_audit,
                  "artifacts": list({a["path"]: a for a in artifacts}.values()),
                  "pocket": {"matched_context": False, "supported_herbicide_positions": []},
                  "blocking_reasons": ["Chemical states lack a validated pH/protonation protocol",
                                       "Experimental controls are homologs, not a same-target Arabidopsis accuracy benchmark",
                                       "Known-control kinetic endpoints do not calibrate predicted pIC50 or establish function",
                                       "Independent S3P-aware pocket corroboration remains unavailable"],
                  "legend": LEGEND}
        if not matched_control_msa:
            report["blocking_reasons"].append("Control MSAs differ; mutant binding comparisons are confounded by alignment differences")
        write(self.root / "predictions.json", {"records": predictions, "legend": LEGEND})
        for filename in ("predictions.json", "calibration_inputs.json", "chemical_state_audit.json", "evidence.json"):
            report["artifacts"].append({"path": filename, "sha256": sha(self.root / filename)})
        from .binding_calibration import build_binding_report
        binding = build_binding_report(self.root, self.root / "binding_dataset.json" if (self.root / "binding_dataset.json").exists() else None,
                                       input_report=report)
        write(self.root / "binding_calibration.json", binding)
        report["binding_calibration"] = {"status": binding["status"], "path": "binding_calibration.json"}
        report["artifacts"].append({"path": "binding_calibration.json", "sha256": sha(self.root / "binding_calibration.json")})
        write(self.root / "calibration_report.json", report)
        rows = ["# EPSPS Calibration", "", "Status: NEEDS_REASSESSMENT. No new candidates advanced.", "", LEGEND, "",
                "| Comparison | Subject | Context | Global RMSD (A) | Pocket RMSD (A) | TM query | Failed checks |",
                "|---|---|---|---:|---:|---:|---|"]
        rows += [f"| {r['kind']} | {r['subject']} | {r['context']} | {r['global_ca_rmsd_angstrom']:.3f} | "
                 f"{format(r['active_site_rmsd_angstrom'], '.3f') if r['active_site_rmsd_angstrom'] is not None else 'NA'} | "
                 f"{r['query_tm_score']:.3f} | {', '.join(r['failed_metrics']) or 'none'} |"
                 for r in comparisons]
        rows += ["", "Table legend: " + LEGEND, "", "## Remaining Requirements", ""] + ["- " + s for s in report["blocking_reasons"]]
        rows += ["", "## Binding-Control Diagnostics", "",
                 "| Ligand | Mean mutant minus WT predicted pIC50 | Predicted weaker in mutant |",
                 "|---|---:|---|"]
        rows += [f"| {r['ligand']} | {r['mean_mutant_minus_wt_pIC50']:.3f} | {r['predicted_weaker_in_mutant']} |" for r in controls]
        rows += ["", "Table legend: means use three seeds per E. coli sequence. A negative difference means weaker predicted binding. "
                 "The G96A control is reported to have glyphosate insensitivity and reduced PEP affinity in its associated study "
                 "(https://www.rcsb.org/structure/1MI4). These direction checks do not equate predicted pIC50 with kinetic constants "
                 "or provide a validated equivalence test. WT and mutant homolog MSA rows matched: " + str(matched_control_msa) + ".",
                 "", "## Contact Audit", "", contact_audit["legend"], "",
                 "Herbicide-only positions under the existing rule: " + str(contact_audit["herbicide_only_positions"]) + ".",
                 "", f"Predictions analysed: {len(predictions)}; new in this run: {manifest.get('new_prediction_count', len(predictions))}; "
                 f"hash-validated reuse: {manifest.get('reused_prediction_count', 0)}."]
        (self.root / "CALIBRATION_REPORT.md").write_text("\n".join(rows) + "\n")
        execution_path = self.root / "execution.json"
        if execution_path.exists():
            state = json.loads(execution_path.read_text())
            complete = {s["stage"] for s in state["stages"] if s["status"] == "COMPLETED"}
            expected = {"msa", *[f"boltz-{seed}" for seed in manifest["seeds"]]}
            if expected <= complete:
                if state.get("failure_type"):
                    state.setdefault("recovered_failures", []).append(state.pop("failure_type"))
                state["status"] = "COMPUTATION_COMPLETE"
                state["analysis"] = {"status": "VALIDATED", "prediction_count": len(predictions),
                                     "comparison_count": len(comparisons), "source_sha256": sha(Path(__file__)),
                                     "report_sha256": sha(self.root / "calibration_report.json")}
                write(execution_path, state)
        return report
