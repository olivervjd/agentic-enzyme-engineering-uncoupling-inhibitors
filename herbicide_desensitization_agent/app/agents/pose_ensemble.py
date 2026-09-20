"""Public DiffDock-L execution with explicit protein-only context limitations."""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

from ..backends.calibration import sha, write
from ..chemistry.contacts import residue_contacts_from_pdb_and_sdf


IMAGE = "rbgcsail/diffdock@sha256:1b7bb3adb332fdc9648a0ec53dec2f790cfbb816d7478d2bd93f1cdea3b269f0"
LEGEND = ("Independent DiffDock-L diagnostic poses on Boltz WT protein coordinates, four poses per receptor. "
          "Contacts use a 5 angstrom heavy-atom cutoff and full Arabidopsis sequence numbering (+76). "
          "S3P is deliberately excluded because this backend cannot represent it as a fixed cosubstrate. "
          "These are not context-matched confirmations of Boltz S3P-containing complexes. Raw pose confidence "
          "is not binding affinity or a probability. No diagnostic pose can override the quality/pocket gate.")


class PoseEnsembleAgent:
    def __init__(self, output, cache, timeout=3600):
        self.root, self.cache, self.timeout = Path(output), Path(cache), timeout

    def run(self, request, structures):
        from Bio.PDB import MMCIFParser, PDBIO, Select
        from Bio.SeqUtils import seq1
        from rdkit import Chem

        if self.root.exists():
            raise ValueError("Use a fresh diagnostic docking directory")
        selected = [s for s in structures if s.scores.get("context") in {"native", "herbicide"}]
        if not selected or len({s.model_id for s in selected}) != len(selected):
            raise ValueError("Unique native/herbicide receptor models are required")
        self.root.mkdir(parents=True)
        (self.root / "inputs").mkdir()
        script = Path(__file__).resolve().parents[2] / "examples" / "run_diffdock_smoke.py"
        shutil.copy2(script, self.root / "run_diffdock.py")
        records = []
        native = next(l for l in request.native_ligands if l.name == "phosphoenolpyruvate")
        for model in selected:
            if not re.fullmatch(r"[A-Za-z0-9_-]+", model.model_id):
                raise ValueError("Unsafe docking record ID")
            parsed = MMCIFParser(QUIET=True).get_structure(model.model_id, model.artifact_path)
            sequence = "".join(seq1(r.resname) for r in parsed[0]["A"] if r.id[0] == " ")
            if sequence != request.target.sequence[76:]:
                raise ValueError("Docking receptor is not the expected WT sequence")
            class ProteinOnly(Select):
                def accept_model(self, model):
                    return model.id == 0

                def accept_chain(self, chain):
                    return chain.id == "A"

                def accept_residue(self, residue):
                    return residue.id[0] == " "
            receptor = self.root / "inputs" / f"{model.model_id}.pdb"
            writer = PDBIO()
            writer.set_structure(parsed)
            writer.save(str(receptor), ProteinOnly())
            ligand = request.herbicide if model.scores["context"] == "herbicide" else native
            if ligand.structure_format != "SMILES" or Chem.MolFromSmiles(ligand.structure) is None:
                raise ValueError("Diagnostic ligand requires valid SMILES")
            records.append({"id": model.model_id, "ligand": ligand.name, "smiles": ligand.structure,
                            "receptor": str(receptor.relative_to(self.root)), "receptor_sha256": sha(receptor),
                            "source_sha256": sha(model.artifact_path), "source": model.artifact_path})
        with (self.root / "inputs/diffdock.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["complex_name", "protein_path", "ligand_description", "protein_sequence"])
            writer.writeheader()
            for row in records:
                writer.writerow({"complex_name": row["id"], "protein_path": "/validation/" + row["receptor"],
                                 "ligand_description": row["smiles"], "protein_sequence": ""})
        container = "epsps-diagnostic-" + sha(self.root / "inputs/diffdock.csv")[:12]
        command = ["docker", "run", "--rm", "--name", container, "--network", "none", "--shm-size=4g",
                   "-e", "PYTHONPATH=/home/appuser/DiffDock", "-e", "OMP_NUM_THREADS=4", "-e", "DIFFDOCK_THREADS=4",
                   "-e", "DIFFDOCK_DEVICE=cpu", "-e", "DIFFDOCK_SEED=101", "-v", f"{self.root.resolve()}:/validation",
                   "-v", f"{self.cache.resolve()}/workdir:/home/appuser/DiffDock/workdir:ro",
                   "-v", f"{self.cache.resolve()}/torch:/home/appuser/.cache/torch:ro", IMAGE,
                   "micromamba", "run", "-n", "diffdock", "python", "-u", "/validation/run_diffdock.py"]
        state = {"status": "RUNNING", "command": command, "records": records, "legend": LEGEND}
        write(self.root / "execution.json", state)
        try:
            with (self.root / "run.log").open("w") as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=self.timeout)
            if result.returncode:
                raise RuntimeError("DiffDock diagnostic failed; see run.log")
            poses = self.collect(records)
            report = {"status": "COMPLETED_DIAGNOSTIC_ONLY", "matched_context": False,
                      "records": poses, "legend": LEGEND}
            write(self.root / "pose_diagnostics.json", report)
            state["status"] = "COMPLETED"
            return report
        except BaseException:
            state["status"] = "FAILED"
            subprocess.run(["docker", "stop", "--time", "5", container], capture_output=True, timeout=30)
            raise
        finally:
            write(self.root / "execution.json", state)

    def collect(self, records):
        from rdkit import Chem

        poses = []
        for row in records:
            receptor = self.root / row["receptor"]
            if sha(receptor) != row["receptor_sha256"] or sha(row["source"]) != row["source_sha256"]:
                raise ValueError("Docking receptor changed during execution")
            paths = list((self.root / "diffdock" / row["id"]).glob("rank*_confidence*.sdf"))
            ranks = []
            for path in paths:
                match = re.fullmatch(r"rank(\d+)_confidence([-+0-9.eE]+)\.sdf", path.name)
                if not match:
                    raise ValueError("Unrecognized pose filename")
                rank, confidence = int(match[1]), float(match[2])
                ranks.append(rank)
                mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=True)), None)
                expected = Chem.MolFromSmiles(row["smiles"])
                if mol is None or Chem.MolToSmiles(mol) != Chem.MolToSmiles(expected):
                    raise ValueError("Docked ligand chemistry mismatch")
                if not math.isfinite(confidence) or any(not math.isfinite(float(x)) for xyz in mol.GetConformer().GetPositions() for x in xyz):
                    raise ValueError("Non-finite docking prediction")
                contacts = residue_contacts_from_pdb_and_sdf(receptor.read_text(), Chem.MolToMolBlock(mol),
                                                              cutoff_angstrom=5.0, protein_chain="A", residue_offset=76)
                poses.append({"model_id": row["id"], "ligand": row["ligand"], "rank": rank,
                              "raw_confidence": confidence, "contacts": contacts,
                              "path": str(path.relative_to(self.root)), "sha256": sha(path)})
            if sorted(ranks) != [1, 2, 3, 4]:
                raise ValueError("Expected four distinct ranked poses per receptor")
        return poses
