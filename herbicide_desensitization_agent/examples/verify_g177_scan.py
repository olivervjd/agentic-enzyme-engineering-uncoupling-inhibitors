"""Read-only scientific artifact integrity audit for the G177 scan."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from Bio.PDB import MMCIFParser, Superimposer
from Bio.SeqUtils import seq1
from rdkit import Chem


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(root):
    manifest = json.loads((root / "manifest.json").read_text())
    results = json.loads((root / "results.json").read_text())
    if not results["complete"] or results["completed"] != 120:
        raise ValueError("Refusing final audit of incomplete campaign")
    inputs = {x["id"]: x for x in manifest["records"]}
    structures, comparisons = {}, []
    common_homologs = None
    for row in results["replicates"]:
        registered = inputs[row["id"]]
        if digest(root / registered["request"]) != registered["request_sha256"]:
            raise ValueError("Request changed")
        msa_path = root / registered["msa"]
        if digest(msa_path) != registered["msa_sha256"]:
            raise ValueError("MSA changed")
        msa = msa_path.read_text().splitlines(keepends=True)
        headers = [i for i, line in enumerate(msa) if line.startswith(">")]
        query = "".join(line.strip() for line in msa[headers[0]+1:headers[1]])
        tail = "".join(msa[headers[1]:])
        if common_homologs is None:
            common_homologs = tail
        if tail != common_homologs or query != registered["mature_sequence"]:
            raise ValueError("MSA background differs or query is incorrect")
        paths = [root / a["path"] for a in row["artifacts"]]
        for path, artifact in zip(paths, row["artifacts"]):
            if digest(path) != artifact["sha256"]:
                raise ValueError("Output artifact changed")
        cif = next(p for p in paths if p.suffix == ".cif")
        model = MMCIFParser(QUIET=True).get_structure(row["id"], str(cif))[0]
        residues = [r for r in model["A"] if r.id[0] == " "]
        if "".join(seq1(r.resname) for r in residues) != registered["mature_sequence"]:
            raise ValueError("Incorrect output sequence")
        if len(residues) != 444 or any("CA" not in r for r in residues):
            raise ValueError("Incomplete protein chain")
        if seq1(residues[100].resname) != row["amino_acid"]:
            raise ValueError("Wrong amino acid at canonical G177")
        chemistry = manifest["chemical_states"][row["ligand"]]
        for chain, smiles in (("B", chemistry["encoded_smiles"]), ("C", chemistry["s3p_smiles"])):
            molecule = Chem.MolFromSmiles(smiles)
            expected = Counter(a.GetSymbol().upper() for a in molecule.GetAtoms() if a.GetAtomicNum() > 1)
            actual = Counter(a.element.upper() for a in model[chain].get_atoms() if a.element.upper() != "H")
            if actual != expected:
                raise ValueError(f"Ligand heavy-atom inventory mismatch: {row['id']} / {chain}")
        structures[(row["variant"], row["ligand"], row["seed"])] = residues
    for (variant, ligand, seed), residues in structures.items():
        wt = structures[("WT", ligand, seed)]
        sup = Superimposer()
        sup.set_atoms([r["CA"] for r in wt], [r["CA"] for r in residues])
        comparisons.append({"variant": variant, "ligand": ligand, "seed": seed,
            "global_ca_rmsd_angstrom_vs_matched_wt": float(sup.rms)})
    return {"status": "PASSED", "validated_jobs": len(structures),
        "checks": ["120 unique complete jobs", "input and output SHA256 integrity",
            "identical homolog MSA background and correct mutant query", "444-residue exact output sequence",
            "correct canonical residue 177 / mature residue 101", "PEP/glyphosate and S3P heavy-atom inventories"],
        "limitations": ["Heavy-atom inventories do not prove protonation or stereochemistry",
            "Global CA RMSD is a structural diagnostic, not proof of function retention"],
        "structure_comparisons": comparisons}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.root), indent=2, allow_nan=False))
