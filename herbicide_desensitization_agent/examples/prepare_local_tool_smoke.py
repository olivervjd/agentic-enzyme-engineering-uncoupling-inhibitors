"""Prepare matched WT EPSPS inputs for upstream DiffDock-L and Boltz-2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def prepare(inputs: Path, reference: Path, output: Path) -> dict:
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.SeqUtils import seq1
    from rdkit import Chem

    sequence = "".join(line.strip() for line in (inputs / "P05466.fasta").read_text().splitlines()
                       if not line.startswith(">"))[76:]
    if len(sequence) != 444:
        raise ValueError("Expected the validated 444-residue mature P05466 sequence")
    ligands = {}
    for name in ("pep", "glyphosate", "s3p"):
        data = json.loads((inputs / f"{name}.json").read_text())
        value = data["PropertyTable"]["Properties"][0]["SMILES"]
        if Chem.MolFromSmiles(value) is None:
            raise ValueError(f"Invalid {name} SMILES")
        ligands[name] = value

    structure = MMCIFParser(QUIET=True).get_structure("epsps", reference)
    residues = [residue for residue in structure[0]["A"] if residue.id[0] == " "]
    if "".join(seq1(residue.resname) for residue in residues) != sequence:
        raise ValueError("Reference protein chain A does not match the expected WT sequence")

    class ProteinOnly(Select):
        def accept_model(self, model):
            return model.id == 0

        def accept_chain(self, chain):
            return chain.id == "A"

        def accept_residue(self, residue):
            return residue.id[0] == " "

    output.mkdir(parents=True, exist_ok=True)
    boltz = output / "boltz"
    boltz.mkdir(exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(output / "epsps.pdb"), ProteinOnly())
    rows = []
    for name in ("pep", "glyphosate"):
        record = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": sequence, "msa": "empty"}},
                {"ligand": {"id": "B", "smiles": ligands[name]}},
                {"ligand": {"id": "C", "smiles": ligands["s3p"]}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
        (boltz / f"epsps-{name}.yaml").write_text(json.dumps(record, indent=2) + "\n")
        rows.append({"complex_name": f"epsps-{name}", "protein_path": "/validation/inputs/epsps.pdb",
                     "ligand_description": ligands[name], "protein_sequence": ""})
    with (output / "diffdock.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "target": "AT2G45300", "reference": "wild-type", "mature_sequence_length": len(sequence),
        "residue_numbering_offset": 76, "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        "ligands": ligands, "msa": "query-only",
        "legend": "Execution smoke test, not a mutation campaign or evidence of retained enzyme function. "
                  "DiffDock-L docks PEP or glyphosate separately to protein-only EPSPS chain A; it excludes S3P. "
                  "Boltz-2 cofolds EPSPS with PEP or glyphosate plus S3P and requests affinity for chain B only. "
                  "These context differences preclude treating the tools as equivalent affinity measurements. "
                  "Both primary ligands use the same stored PubChem SMILES without pH-specific protonation enumeration.",
    }
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.inputs, args.reference, args.output), indent=2))


if __name__ == "__main__":
    main()
