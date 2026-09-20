"""Enumerate EPSPS ligand microstate hypotheses at the published WT assay pH.

No population or bound-state validation is inferred from enumeration. The chosen
campaign states are explicit, reversible modeling assumptions; all alternatives
remain in the audit for sensitivity analysis.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path


CAMPAIGN_STATES = {
    "pep": "C=C(OP(=O)([O-])[O-])C(=O)[O-]",
    "s3p": "O=C([O-])C1=C[C@@H](OP(=O)([O-])[O-])[C@@H](O)[C@H](O)C1",
    "glyphosate": "O=C([O-])C[NH2+]CP(=O)([O-])[O-]",
}


def enumerate_states(ligands):
    from dimorphite_dl import protonate_smiles
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    if set(ligands) != set(CAMPAIGN_STATES):
        raise ValueError("Exactly PEP, S3P and glyphosate are required")
    uncharger = rdMolStandardize.Uncharger()
    def neutral_identity(mol):
        return Chem.MolToSmiles(uncharger.uncharge(mol), isomericSmiles=True)
    result = {}
    for name, smiles in ligands.items():
        original = Chem.MolFromSmiles(smiles)
        if original is None:
            raise ValueError("Invalid input molecule")
        raw = protonate_smiles(smiles, ph_min=7.4, ph_max=7.4, precision=1.0, max_variants=128)
        variants = []
        for state in sorted(set(raw)):
            mol = Chem.MolFromSmiles(state)
            if mol is None or neutral_identity(mol) != neutral_identity(original):
                raise ValueError("Protonation enumeration changed covalent structure or stereochemistry")
            canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
            variants.append({"state_id": name + ":" + hashlib.sha256(canonical.encode()).hexdigest()[:12],
                             "smiles": canonical, "formal_charge": Chem.GetFormalCharge(mol),
                             "formula": rdMolDescriptors.CalcMolFormula(mol),
                             "inchi_key": Chem.MolToInchiKey(mol),
                             "population": None, "covalent_identity_and_stereochemistry_retained": True})
        selected_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(CAMPAIGN_STATES[name]), isomericSmiles=True)
        selected = [v for v in variants if v["smiles"] == selected_smiles]
        if len(selected) != 1:
            raise ValueError("Preregistered campaign hypothesis not returned by enumerator")
        result[name] = {"input_smiles": smiles, "microstates": variants, "selected_campaign_state": selected[0],
                        "selection_basis": "Explicit carboxylate and doubly deprotonated phosphate/phosphonate hypothesis; glyphosate amine protonated",
                        "selection_status": "MODELING_ASSUMPTION_NOT_VALIDATED_MICROSTATE"}
    return {"schema_version": "1.0.0", "target_agi": "AT2G45300", "assay_pH": 7.4,
            "assay_source": "https://doi.org/10.1016/j.csbj.2022.03.020", "assay_source_section": "4.4",
            "enumeration": {"model": "Dimorphite-DL", "version": version("dimorphite_dl"),
                            "rdkit_version": rdBase.rdkitVersion, "ph_min": 7.4, "ph_max": 7.4,
                            "precision": 1.0, "max_variants": 128,
                            "source": "https://doi.org/10.1186/s13321-019-0336-9"},
            "ligands": result,
            "protonation_experimentally_validated": False,
            "limitations": ["Group-based empirical enumeration, not molecule-specific measured pKa or bound-state validation",
                            "Do not choose the best-scoring microstate after candidate inspection",
                            "Run sensitivity analysis on retained alternatives when conclusions depend on charge",
                            "Protein protonation, catalytic waters and ionic strength need separate preparation",
                            "Downstream tools may standardize ionization; compare actual encoded states, not requested SMILES",
                            "Assay pH 7.4 does not imply crystallization pH 7.5 or historical predictions had matching chemistry"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Calibration/input manifest containing ligands")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source = args.input.read_bytes()
    report = enumerate_states(json.loads(source)["ligands"])
    report["input_manifest_sha256"] = hashlib.sha256(source).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({name: {"states": len(r["microstates"]), "selected": r["selected_campaign_state"]}
                      for name, r in report["ligands"].items()}))


if __name__ == "__main__":
    main()
