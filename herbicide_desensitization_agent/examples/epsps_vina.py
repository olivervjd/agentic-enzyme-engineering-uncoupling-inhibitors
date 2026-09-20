"""Replicated Vina pose diagnostics with explicitly retained EPSPS co-substrates.

Vina scores are empirical scoring-function values in kcal/mol, not measured
binding constants or thermodynamic free energies. Preparation fails closed when
atom mapping, chemical graph, or required co-substrate coordinates are missing.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata
import io
import json
import math
import time
from pathlib import Path

SEEDS = (211, 223, 227)
SUPPORTED_CONTEXTS = {"herbicide": {"glyphosate", "s3p"}, "native": {"pep", "s3p"},
                      "native_s3p": {"pep", "s3p"}, "s3p_only": {"s3p"}}
LIMITATIONS = [
    "Vina scores are empirical pose-scoring values, not thermodynamic delta-G or calibrated affinity.",
    "Protein side-chain protonation uses Meeko residue templates, not a validated pH-dependent ensemble.",
    "Non-ligand waters and buffer components are omitted and enumerated; catalytic water requirements are not validated.",
    "Required co-substrate heavy atoms are retained as a rigid receptor component; flexibility and polarization are not modeled.",
    "Vina's empirical scoring function does not explicitly model Coulomb electrostatics; highly charged phosphate/phosphonate scores need independent calibration.",
    "The box is centred on a recorded ligand site. This is targeted redocking/cross-method docking, not blind binding-site discovery.",
    "Same ligand graph/state is preserved, but preparation assumptions differ from complex prediction; pose agreement is not direct binding validation.",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def versions(require_vina=False):
    names = ["rdkit", "meeko", "gemmi", "biopython", "numpy", "scipy"] + (["vina"] if require_vina else [])
    return {name: importlib.metadata.version(name) for name in names}


def heavy_atoms(residue):
    result = {}
    for item in residue:
        atom = item.selected_child if item.is_disordered() == 2 else item
        if atom.element.upper() in {"H", "D"}:
            continue
        if atom.name in result:
            raise ValueError("Duplicate ligand atom name")
        xyz = tuple(float(x) for x in atom.coord)
        if len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
            raise ValueError("Nonfinite ligand coordinates")
        result[atom.name] = {"element": atom.element.upper(), "xyz": xyz}
    return result


def boltz_atom_names(smiles):
    """Exact structure-only Boltz 2.2.1 naming, without affinity standardization."""
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid input SMILES")
    molecule = Chem.AddHs(molecule)
    ranks = Chem.CanonicalRankAtoms(molecule)
    Chem.AssignStereochemistry(molecule, force=True, cleanIt=True)
    return {atom.GetIdx(): atom.GetSymbol().upper() + str(rank + 1)
            for atom, rank in zip(molecule.GetAtoms(), ranks) if atom.GetAtomicNum() != 1}


def molecule_with_coordinates(smiles, atoms, atom_map):
    """Attach validated heavy-atom coordinates to the declared chemical graph."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Expected a single valid ligand graph")
    if set(atom_map) != set(range(mol.GetNumAtoms())) or set(atom_map.values()) != set(atoms):
        raise ValueError("Ligand atom mapping must cover each and every heavy atom exactly once")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for atom in mol.GetAtoms():
        name = atom_map[atom.GetIdx()]
        if atom.GetSymbol().upper() != atoms[name]["element"]:
            raise ValueError("Mapped ligand element mismatch")
        conformer.SetAtomPosition(atom.GetIdx(), atoms[name]["xyz"])
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4)); info.SetResidueName("LIG"); info.SetResidueNumber(1)
        info.SetChainId("L"); info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    mol.AddConformer(conformer)
    for bond in mol.GetBonds():
        a, b = conformer.GetAtomPosition(bond.GetBeginAtomIdx()), conformer.GetAtomPosition(bond.GetEndAtomIdx())
        distance = (a-b).Length()
        if not .7 <= distance <= 2.3:
            raise ValueError("Mapped chemical bond has an implausible heavy-atom distance")
    expected_stereo = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    from_coords = Chem.Mol(mol)
    Chem.RemoveStereochemistry(from_coords)
    Chem.AssignAtomChiralTagsFromStructure(from_coords, replaceExistingTags=True)
    Chem.AssignStereochemistry(from_coords, cleanIt=True, force=True)
    observed_stereo = dict(Chem.FindMolChiralCenters(from_coords, includeUnassigned=True))
    for idx, desired in expected_stereo.items():
        if desired != "?" and observed_stereo.get(idx) != desired:
            raise ValueError("Ligand coordinate stereochemistry disagrees with declared chemical state")
    return mol


def ccd_atom_map(structure_path, component, smiles, atoms):
    """Map the deposited CCD graph to a declared protonation hypothesis.

    Protonation/charges remain those in SMILES. Element and bond-order graph
    equivalence plus 3D stereochemistry are required; atom order is never used.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from rdkit import Chem
    cif = MMCIF2Dict(str(structure_path))
    required = ["_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
                "_chem_comp_bond.comp_id", "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2", "_chem_comp_bond.value_order"]
    if any(k not in cif for k in required):
        raise ValueError("Experimental atom mapping requires deposited CCD atom/bond tables or an explicit validated atom_map")
    ids, index = [], {}
    graph = Chem.RWMol()
    for comp, name, element in zip(*(cif[k] for k in required[:3])):
        if comp == component and element.upper() not in {"H", "D"}:
            index[name] = graph.AddAtom(Chem.Atom(element.title()))
            ids.append(name)
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for comp, a, b, order in zip(*(cif[k] for k in required[3:])):
        if comp == component and a in index and b in index:
            if order.upper() not in orders:
                raise ValueError("Unsupported CCD bond order")
            graph.AddBond(index[a], index[b], orders[order.upper()])
    target = graph.GetMol()
    target.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(target)
    template = Chem.MolFromSmiles(smiles)
    if template is None or template.GetNumAtoms() != len(ids) or set(ids) != set(atoms):
        raise ValueError("CCD, deposited pose and declared ligand have different heavy-atom sets")
    query = Chem.RWMol(template)
    for atom in template.GetAtoms():
        query.ReplaceAtom(atom.GetIdx(), Chem.AtomFromSmarts(f"[#{atom.GetAtomicNum()}]"))
    matches = target.GetSubstructMatches(query.GetMol(), uniquify=False, useChirality=False, maxMatches=10000)
    mappings = []
    for match in matches:
        mapping = {i: ids[j] for i, j in enumerate(match)}
        try:
            molecule_with_coordinates(smiles, atoms, mapping)
        except ValueError:
            continue
        mappings.append(mapping)
    if not mappings:
        raise ValueError("No chemically and stereochemically consistent experimental ligand atom mapping")
    mappings.sort(key=lambda m: tuple(m[i] for i in sorted(m)))
    return mappings[0], len(mappings)


def select_ligand(model, spec):
    chain_id = spec.get("chain", spec.get("id"))
    choices = []
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] == " " or residue.resname in {"HOH", "WAT"}:
                continue
            if spec.get("component") and residue.resname != spec["component"]:
                continue
            if spec.get("residue_number") is not None and residue.id[1] != int(spec["residue_number"]):
                continue
            choices.append(residue)
    if len(choices) != 1:
        raise ValueError(f"Ligand {spec['name']} selection is absent or ambiguous ({len(choices)} residues)")
    return choices[0]


def validate_context(record, target_name):
    if record.get("mode") not in {"structure", "experimental"}:
        raise ValueError("Docking accepts structure-only or experimental coordinates, never affinity-standardized complexes")
    required = SUPPORTED_CONTEXTS.get(record.get("context"))
    if required is None:
        raise ValueError("Explicit supported ligand-bound context is required; apo docking requires a separately registered box/control protocol")
    names = [x["name"] for x in record.get("ligands", [])]
    if len(names) != len(set(names)) or set(names) != required or target_name not in names:
        raise ValueError("All required context ligands must be explicitly represented exactly once")
    for spec in record["ligands"]:
        if spec.get("standardized_by_affinity_parser"):
            raise ValueError("Affinity-standardized ligand is not matched structure-only input")
        if not (spec.get("model_smiles") or spec.get("smiles")):
            raise ValueError("Explicit ligand SMILES required")
    return required - {target_name}


def protein_pdb(model, chains, expected_sequence=None):
    from Bio.PDB import PDBIO, Structure, Model, Chain
    from Bio.PDB.Polypeptide import is_aa
    from Bio.SeqUtils import seq1
    protein = Structure.Structure("receptor")
    clean = Model.Model(0)
    protein.add(clean)
    sequence = ""
    source_residues = []
    for chain_id in chains:
        if chain_id not in model:
            raise ValueError("Required protein chain missing")
        target_chain = Chain.Chain(chain_id)
        clean.add(target_chain)
        for residue in model[chain_id]:
            if not is_aa(residue, standard=False):
                continue
            if residue.is_disordered() == 2 or not is_aa(residue, standard=True):
                raise ValueError("Ambiguous or modified protein residue requires an explicit preparation policy")
            cloned = copy.deepcopy(residue)
            for atom in list(cloned):
                cloned.detach_child(atom.id)
            for item in residue:
                selected = item.selected_child if item.is_disordered() == 2 else item
                if selected.element.upper() in {"H", "D"}:
                    continue
                atom = selected.copy()
                atom.set_altloc(" ")
                atom.disordered_flag = 0
                cloned.add(atom)
            target_chain.add(cloned)
            sequence += seq1(residue.resname)
            source_residues.append(f"{chain_id}:{residue.id[1]}{residue.id[2].strip()}")
    if not source_residues or (expected_sequence is not None and sequence != expected_sequence):
        raise ValueError("Protein sequence missing or differs from registered input")
    writer = PDBIO(); writer.set_structure(protein)
    stream = io.StringIO(); writer.save(stream)
    return stream.getvalue(), source_residues, sequence


def prepare_ligand(molecule, chain, residue_name="LIG"):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    with_hydrogens = Chem.AddHs(molecule, addCoords=True)
    for atom in with_hydrogens.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            info = Chem.AtomPDBResidueInfo()
            info.SetName(("H"+str(atom.GetIdx()+1)).rjust(4))
        info.SetChainId(chain); info.SetResidueName(residue_name); info.SetResidueNumber(1)
        atom.SetMonomerInfo(info)
    setups = MoleculePreparation(charge_model="gasteiger").prepare(with_hydrogens)
    if len(setups) != 1:
        raise ValueError("Ligand preparation yielded an ambiguous setup ensemble")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        raise ValueError("Ligand PDBQT preparation failed: " + str(error))
    if any(not math.isfinite(float(atom.charge)) for atom in setups[0].atoms if not atom.is_ignore):
        raise ValueError("Nonfinite prepared atomic charge")
    # Atom names must survive preparation, including their coordinates.
    expected = {a.GetPDBResidueInfo().GetName().strip(): molecule.GetConformer().GetAtomPosition(a.GetIdx())
                for a in molecule.GetAtoms()}
    seen = {line[12:16].strip(): tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))
            for line in pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() in expected}
    if set(seen) != set(expected):
        raise ValueError("Preparation dropped or renamed required ligand heavy atoms")
    if any(math.dist(seen[name], tuple(expected[name])) > .002 for name in expected):
        raise ValueError("Preparation moved ligand heavy-atom coordinates")
    return pdbqt


def validate_protein_coordinates(source, prepared, tolerance=.002):
    """Permit only an exact terminal carboxylate oxygen name permutation."""
    if set(source) - set(prepared):
        raise ValueError("Protein preparation omitted source heavy atoms")
    remapped=[]
    for key,xyz in source.items():
        if math.dist(xyz,prepared[key])<=tolerance:continue
        other=(*key[:2],"OXT" if key[2]=="O" else "O")
        if key[2] in {"O","OXT"} and other in source and other in prepared and math.dist(xyz,prepared[other])<=tolerance and math.dist(source[other],prepared[key])<=tolerance:
            remapped.append({"chain":key[0],"residue":key[1].strip(),"source_atom":key[2],"prepared_atom":other[2],
                             "reason":"exact terminal carboxylate oxygen name permutation; coordinates retained"})
        else:
            raise ValueError("Protein preparation moved source heavy atoms")
    return remapped


def combine_rigid_receptor(protein_pdbqt, cosubstrate_pdbqts):
    lines = [line for line in protein_pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    protein_count = len(lines)
    counts = {}
    for name, pdbqt in cosubstrate_pdbqts.items():
        added = [line for line in pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
        if not added:
            raise ValueError("Required rigid co-substrate has no prepared atoms")
        counts[name] = len(added)
        lines.extend(added)
    if protein_count == 0:
        raise ValueError("Prepared protein receptor has no atoms")
    result = "\n".join(line[:6]+f"{i:5d}"+line[11:] for i, line in enumerate(lines,1)) + "\n"
    if any(line.startswith(("ROOT", "BRANCH", "TORSDOF")) for line in result.splitlines()):
        raise ValueError("Rigid co-substrate cannot contain ligand torsion records")
    return result, {"protein_atom_count":protein_count,"rigid_cosubstrate_atom_counts":counts,"total_atom_count":len(lines)}


def prepare_record(record, campaign_dir, output, target_name, expected_sequence=None):
    from Bio.PDB import MMCIFParser
    from meeko import Polymer, MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    required_cosubstrates = validate_context(record, target_name)
    path = (Path(campaign_dir)/record["structure"]).resolve()
    if record.get("structure_sha256") and digest(path) != record["structure_sha256"]:
        raise ValueError("Input coordinate artifact hash mismatch")
    model = MMCIFParser(QUIET=True).get_structure(record["id"], path)[0]
    molecules, ligand_records = {}, []
    ligand_residues = set()
    for spec in record["ligands"]:
        residue = select_ligand(model, spec)
        ligand_residues.add((residue.parent.id,residue.id))
        atoms = heavy_atoms(residue)
        smiles = spec.get("model_smiles") or spec["smiles"]
        if spec.get("atom_map"):
            mapping = {int(k):v for k,v in spec["atom_map"].items()}; count=1; mapping_method="explicit atom map with graph/geometry validation"
        elif record["mode"] == "structure":
            mapping = boltz_atom_names(smiles); count=1; mapping_method="Boltz 2.2.1 AddHs+CanonicalRankAtoms"
        else:
            mapping,count=ccd_atom_map(path,spec.get("component",residue.resname),smiles,atoms)
            mapping_method="deposited CCD graph isomorphism and coordinate stereochemistry"
        mol = molecule_with_coordinates(smiles,atoms,mapping)
        molecules[spec["name"]] = mol
        writer=Chem.SDWriter(str(output/f"{spec['name']}.sdf")); writer.write(mol); writer.close()
        ligand_records.append({"name":spec["name"],"source_chain":residue.parent.id,"source_residue":residue.id[1],
             "component":residue.resname,"smiles":smiles,"canonical_smiles":Chem.MolToSmiles(mol),
             "formal_charge":Chem.GetFormalCharge(mol),"atom_map":mapping,"mapping_method":mapping_method,
             "equivalent_graph_mappings":count,"heavy_atom_count":mol.GetNumAtoms(),
             "stereochemistry":Chem.FindMolChiralCenters(mol,includeUnassigned=True),
             "role":"docked" if spec["name"]==target_name else "rigid_cosubstrate",
             "prepared_pdbqt_chain":"L" if spec["name"]==target_name else "Q", "prepared_pdbqt_residue":1})
    chains = record.get("protein_chains", ["A"])
    pdb, source_residues, observed_sequence = protein_pdb(model,chains,expected_sequence)
    (output/"protein.pdb").write_text(pdb)
    polymer = Polymer.from_pdb_string(pdb,mk_prep=MoleculePreparation(charge_model="gasteiger"),allow_bad_res=False)
    valid=set(polymer.get_valid_monomers())
    if valid != set(source_residues):
        raise ValueError("Protein preparation omitted registered residues")
    receptor, flexible = PDBQTWriterLegacy.write_from_polymer(polymer)
    protein_templates = {key: monomer.residue_template_key for key,monomer in polymer.get_valid_monomers().items()}
    source_heavy = {(line[21],line[22:27],line[12:16].strip()): tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))
                    for line in pdb.splitlines() if line.startswith(("ATOM  ","HETATM"))}
    prepared_heavy = {(line[21],line[22:27],line[12:16].strip()): tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))
                      for line in receptor.splitlines() if line.startswith(("ATOM  ","HETATM")) and not line.split()[-1].startswith("H")}
    coordinate_remapping=validate_protein_coordinates(source_heavy,prepared_heavy)
    added_heavy = [{"chain":key[0],"residue":key[1].strip(),"atom":key[2],"xyz":prepared_heavy[key]}
                   for key in sorted(set(prepared_heavy)-set(source_heavy))]
    if flexible:
        raise ValueError("Unregistered flexible receptor residues")
    prepared={name:prepare_ligand(mol,"L" if name==target_name else "Q") for name,mol in molecules.items()}
    if not required_cosubstrates.issubset(prepared):
        raise ValueError("Required co-substrate was lost during preparation")
    rigid, counts=combine_rigid_receptor(receptor,{name:prepared[name] for name in required_cosubstrates})
    (output/"receptor.pdbqt").write_text(rigid); (output/"ligand.pdbqt").write_text(prepared[target_name])
    target=molecules[target_name]
    coordinates=[tuple(target.GetConformer().GetAtomPosition(i)) for i in range(target.GetNumAtoms())]
    lo=[min(p[i] for p in coordinates) for i in range(3)]; hi=[max(p[i] for p in coordinates) for i in range(3)]
    center=[(a+b)/2 for a,b in zip(lo,hi)]; box=[max(20.0,b-a+12.0) for a,b in zip(lo,hi)]
    excluded=[]
    for chain in model:
        for residue in chain:
            if (chain.id,residue.id) in ligand_residues or (chain.id in chains and residue.id[0]==" "):
                continue
            excluded.append({"chain":chain.id,"residue":residue.id[1],"component":residue.resname,
                             "reason":"not part of declared protein or ligand context"})
    audit={"record_id":record["id"],"subject":record["subject"],"context":record["context"],"docked_ligand":target_name,
           "mode":"experimental_redocking" if record["mode"]=="experimental" else "cross_method_predicted_receptor_docking",
           "source_structure":str(path),"source_sha256":digest(path),"ligands":ligand_records,
           "protein_chains":chains,"protein_sequence":observed_sequence,"protein_residues":source_residues,
           "canonical_offset":record.get("canonical_offset"),"numbering_status":"declared campaign offset" if record.get("canonical_offset") is not None else "canonical numbering unvalidated; structure numbering retained",
           "protein_protonation":"Meeko default residue templates; not experimentally pH-validated",
           "protein_residue_templates":protein_templates,"protein_heavy_atoms_added_by_templates":added_heavy,
           "protein_atom_name_permutations":coordinate_remapping,
           "alternate_location_policy":"BioPDB highest-occupancy selected conformer per atom; alternate-conformation ensemble not sampled",
           "charge_model":"Gasteiger","rigid_cosubstrates":sorted(required_cosubstrates),**counts,
           "box_center_angstrom":center,"box_size_angstrom":box,"excluded_components":excluded,
           "chemical_equivalence":"Ligand graphs/states matched; protein/water assumptions non-equivalent",
           "versions":versions(),"limitations":LIMITATIONS,
           "files":[{"path":p.name,"sha256":digest(p)} for p in sorted(output.iterdir()) if p.is_file()]}
    save(output/"preparation.json",audit)
    return audit,target


def pose_molecule(reference, block):
    from rdkit import Chem
    molecule=Chem.Mol(reference)
    names={a.GetPDBResidueInfo().GetName().strip():a.GetIdx() for a in molecule.GetAtoms()}
    seen=set(); conformer=Chem.Conformer(molecule.GetNumAtoms())
    for line in block.splitlines():
        if not line.startswith(("ATOM  ","HETATM")):
            continue
        name=line[12:16].strip()
        if name not in names:
            continue
        if name in seen:
            raise ValueError("Duplicate atom in docked pose")
        seen.add(name)
        conformer.SetAtomPosition(names[name],tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))))
    if seen != set(names):
        raise ValueError("Docked pose has incomplete atom mapping")
    molecule.RemoveAllConformers();molecule.AddConformer(conformer)
    return molecule


def pose_blocks(pdbqt):
    blocks=[];current=[]
    for line in pdbqt.splitlines():
        if line.startswith("MODEL"):
            current=[]
        elif line.startswith("ENDMDL"):
            blocks.append("\n".join(current));current=[]
        else:
            current.append(line)
    if current and any(line.startswith(("ATOM  ","HETATM")) for line in current):
        blocks.append("\n".join(current))
    return blocks


def cluster_poses(molecules, threshold=2.0):
    from rdkit.Chem import rdMolAlign
    clusters=[]
    assignments=[]
    for i,mol in enumerate(molecules):
        distances=[float(rdMolAlign.CalcRMS(mol,molecules[c[0]],maxMatches=10000)) for c in clusters]
        matches=[(v,j) for j,v in enumerate(distances) if v<=threshold]
        if matches:
            _,selected=min(matches);clusters[selected].append(i)
        else:
            selected=len(clusters);clusters.append([i])
        assignments.append(selected)
    return assignments,clusters


def protein_atom_records(source_path, chains, offset=None):
    """Coordinates remain in the same receptor frame as Vina's returned poses."""
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa
    model=MMCIFParser(QUIET=True).get_structure("contacts",source_path)[0]
    records=[]
    for chain_id in chains:
        for residue in model[chain_id]:
            if not is_aa(residue,standard=True):continue
            atoms=heavy_atoms(residue)
            records.append({"chain":chain_id,"structure_residue":residue.id[1],
                "insertion_code":residue.id[2].strip(),"residue":residue.id[1]+offset if offset is not None else None,
                "amino_acid":residue.resname,"atoms":atoms})
    return records


def pose_contacts(molecule, pose_record, protein, cutoff=5.0):
    import numpy as np
    ligand=[(a.GetPDBResidueInfo().GetName().strip(),tuple(molecule.GetConformer().GetAtomPosition(a.GetIdx())))
            for a in molecule.GetAtoms() if a.GetAtomicNum()!=1]
    ligand_xyz=np.asarray([xyz for _,xyz in ligand],dtype=float)
    rows=[]
    for residue in protein:
        names=list(residue["atoms"])
        protein_xyz=np.asarray([residue["atoms"][name]["xyz"] for name in names],dtype=float)
        distances=np.linalg.norm(protein_xyz[:,None,:]-ligand_xyz[None,:,:],axis=2)
        pairs=[]
        for i,j in zip(*np.where(distances<=cutoff)):
            pairs.append({"protein_atom":names[i],"ligand_atom":ligand[j][0],"distance_angstrom":float(distances[i,j]),
                          "protein_xyz":list(protein_xyz[i]),"ligand_xyz":list(ligand_xyz[j]),
                          "interaction_type":"general proximity contact"})
        if pairs:
            nearest=min(pairs,key=lambda p:p["distance_angstrom"])
            rows.append({**{k:v for k,v in residue.items() if k!="atoms"},
                "variant":pose_record["subject"],"context":pose_record["context"],"ligand":pose_record["ligand"],
                "source_record":pose_record["record_id"],"pose_id":pose_record["pose_id"],
                "seed":pose_record["seed"],"pose_rank":pose_record["pose_rank"],
                "min_distance":nearest["distance_angstrom"],"distance_unit":"angstrom",
                "protein_atom":nearest["protein_atom"],"ligand_atom":nearest["ligand_atom"],
                "interaction_types":["general proximity contact"],"atom_pairs":pairs})
    return rows


def contact_frequencies(contacts, poses):
    from collections import defaultdict
    grouped=defaultdict(list)
    seeds={p["seed"] for p in poses}
    for r in contacts:grouped[(r["chain"],r["structure_residue"],r["insertion_code"],r["ligand"])].append(r)
    frequencies=[]
    for _,records in sorted(grouped.items()):
        first=records[0]
        support={r["seed"] for r in records}
        pose_ids={r["pose_id"] for r in records}
        top_support={r["seed"] for r in records if r["pose_rank"]==1}
        frequencies.append({k:first[k] for k in ("variant","context","ligand","chain","structure_residue","residue","insertion_code","amino_acid")}|{
            "source_record":first["source_record"],"pose_fraction":len(pose_ids)/len(poses),
            "independent_seed_frequency":len(support)/len(seeds),"top_pose_seed_frequency":len(top_support)/len(seeds),
            "contact_frequency":len(top_support)/len(seeds),"frequency_definition":"one vote from top-ranked pose per independent docking seed",
            "seeds_supporting_any_pose":sorted(support),"seeds_supporting_top_pose":sorted(top_support),
            "poses_supporting":sorted(pose_ids),"replicate_count":len(seeds),"pose_count":len(poses),
            "min_distance":min(r["min_distance"] for r in records),"distance_unit":"angstrom",
            "interaction_types":["general proximity contact"],"method":"AutoDock Vina",
            "uncertainty":"Observed sampling fractions; poses within one seed are correlated; no calibrated probability or interval"})
    return frequencies


def dock_prepared(audit, reference, output, seeds=SEEDS, exhaustiveness=8, n_poses=5):
    from vina import Vina
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from ..app.orchestrator.execution import ExecutionJournal
    output=Path(output)
    journal=ExecutionJournal(output/"docking_journal.json")
    rows=[];molecules=[];all_contacts=[]
    protein=protein_atom_records(audit["source_structure"],audit["protein_chains"],audit.get("canonical_offset"))
    try:
        for seed in seeds:
            def run_seed():
                started=time.monotonic()
                vina=Vina(sf_name="vina",cpu=1,seed=int(seed),verbosity=0)
                vina.set_receptor(str(output/"receptor.pdbqt"))
                vina.set_ligand_from_file(str(output/"ligand.pdbqt"))
                vina.compute_vina_maps(center=audit["box_center_angstrom"],box_size=audit["box_size_angstrom"])
                vina.dock(exhaustiveness=exhaustiveness,n_poses=n_poses)
                scores=vina.energies(n_poses=n_poses,energy_range=100.0)
                pdbqt=vina.poses(n_poses=n_poses,energy_range=100.0)
                blocks=pose_blocks(pdbqt)
                if not blocks or len(blocks)!=len(scores):
                    raise ValueError("Vina pose/score counts disagree or no pose returned")
                (output/f"poses-{seed}.pdbqt").write_text(pdbqt)
                sdf=Chem.SDWriter(str(output/f"poses-{seed}.sdf"))
                seed_rows=[]
                for rank,(block,score) in enumerate(zip(blocks,scores),1):
                    mol=pose_molecule(reference,block)
                    value=float(score[0])
                    if not math.isfinite(value):raise ValueError("Nonfinite Vina score")
                    rmsd=float(rdMolAlign.CalcRMS(mol,reference,maxMatches=10000))
                    row={"record_id":audit["record_id"],"subject":audit["subject"],"context":audit["context"],
                         "ligand":audit["docked_ligand"],"seed":seed,"pose_rank":rank,"pose_id":f"{seed}:{rank}",
                         "vina_score_kcal_mol":value,"score_quantity":"Vina empirical score, not thermodynamic delta-G",
                         "reference_pose_rmsd_angstrom":rmsd,"reference_scope":audit["mode"],
                         "rmsd_method":"RDKit symmetry-aware heavy-atom CalcRMS in receptor frame; no ligand alignment",
                         "rigid_cosubstrates":audit["rigid_cosubstrates"],"independent_method":"AutoDock Vina"}
                    mol.SetProp("pose_id",row["pose_id"]);mol.SetDoubleProp("vina_score_kcal_mol",value)
                    contacts=pose_contacts(mol,row,protein)
                    row["contact_residue_count"]=len(contacts)
                    all_contacts.extend(contacts)
                    sdf.write(mol);molecules.append(mol);rows.append(row);seed_rows.append(row)
                sdf.close()
                return {"seed":seed,"poses":seed_rows,"runtime_seconds":time.monotonic()-started}
            journal.execute("docking",run_seed,inputs={"seed":seed,"exhaustiveness":exhaustiveness,"n_poses":n_poses,
                "cpu":1,"scoring_function":"vina","box_center":audit["box_center_angstrom"],"box_size":audit["box_size_angstrom"]},
                input_artifacts=[output/"preparation.json",output/"receptor.pdbqt",output/"ligand.pdbqt"],
                output_artifacts=[output/f"poses-{seed}.pdbqt",output/f"poses-{seed}.sdf"],evidence_status="AVAILABLE")
    finally:
        journal.finalize()
    assignments,clusters=cluster_poses(molecules)
    for row,assignment in zip(rows,assignments):row["pose_cluster"]=assignment
    cluster_summary=[{"cluster":i,"pose_ids":[rows[j]["pose_id"] for j in indexes],
                      "independent_seeds":sorted({rows[j]["seed"] for j in indexes}),
                      "fraction_seeds":len({rows[j]["seed"] for j in indexes})/len(seeds)} for i,indexes in enumerate(clusters)]
    frequencies=contact_frequencies(all_contacts,rows)
    save(output/"pose_contacts.json",all_contacts)
    save(output/"contact_frequencies.json",frequencies)
    result={"status":"COMPLETED","evidence_status":"AVAILABLE","binding_evidence_state":"Insufficient evidence",
            "rows":rows,"pose_contacts":all_contacts,"residue_contact_frequencies":frequencies,
            "contact_definition":"minimum protein–ligand heavy-atom distance <=5.0 angstrom; general proximity only",
            "pose_clusters":cluster_summary,"clustering_threshold_angstrom":2.0,
            "clustering_method":"Greedy seed/rank-order leader clustering; symmetry-aware heavy-atom RMSD, no fitting",
            "preparation":audit,"versions":versions(require_vina=True),"limitations":LIMITATIONS}
    save(output/"docking_results.json",result)
    return result


def jobs_from_manifest(manifest, record_ids=None):
    selected=[]
    for record in manifest["records"]:
        if record.get("mode") not in {"structure","experimental"} or not record.get("structure"):
            continue
        if record_ids and record["id"] not in record_ids:
            continue
        if record.get("context") not in SUPPORTED_CONTEXTS:
            continue
        # Native complexes evaluate PEP with S3P and S3P with PEP independently.
        targets=["pep","s3p"] if record["context"] in {"native","native_s3p"} else ["glyphosate"] if record["context"]=="herbicide" else ["s3p"]
        selected.extend((record,target) for target in targets)
    if not selected:
        raise ValueError("No completed, registered structure-only or experimental ligand-bound records available")
    return selected


def run(campaign_dir, output, *, manifest_path=None, record_ids=None, prepare_only=False):
    from ..app.orchestrator.execution import ExecutionJournal
    campaign_dir=Path(campaign_dir).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError("Fresh docking output directory required; existing computations are never reused")
    output.mkdir(parents=True)
    path=Path(manifest_path) if manifest_path else campaign_dir/"effective_model_inputs.json"
    manifest=json.loads(path.read_text())
    jobs=jobs_from_manifest(manifest,record_ids)
    journal=ExecutionJournal(output/"execution_journal.json")
    records=[]
    for record,target in jobs:
        folder=output/f"{record['id']}--{target}"
        subject=manifest.get("subjects",{}).get(record["subject"],{})
        expected=subject.get("sequence")
        record={**record,"canonical_offset":record.get("canonical_offset",subject.get("offset"))}
        try:
            event=journal.start("docking",inputs={"operation":"receptor_ligand_preparation","record":record["id"],"ligand":target},
                input_artifacts=[path,campaign_dir/record["structure"]])
            try:
                audit,reference=prepare_record(record,campaign_dir,folder,target,expected)
                journal.finish(event,outputs=audit,output_artifacts=list(folder.glob("*")),evidence_status="AVAILABLE")
            except Exception as error:
                journal.fail(event,error)
                raise
            if prepare_only:
                records.append({"id":record["id"],"ligand":target,"execution_status":"PREPARED_ONLY","preparation":str(folder/"preparation.json")})
            else:
                result=journal.execute("docking",lambda:dock_prepared(audit,reference,folder),
                    inputs={"operation":"vina_ensemble","record":record["id"],"ligand":target,
                            "seeds":list(SEEDS),"exhaustiveness":8,"n_poses":5},
                    input_artifacts=[folder/"preparation.json",folder/"receptor.pdbqt",folder/"ligand.pdbqt"],
                    output_artifacts=[folder/"docking_results.json",folder/"docking_journal.json"],
                    evidence_status="AVAILABLE")
                records.append({"id":record["id"],"ligand":target,"execution_status":"COMPLETED","result":str(folder/"docking_results.json"),"poses":len(result["rows"])})
        except Exception as error:
            records.append({"id":record["id"],"ligand":target,"execution_status":"FAILED","error_type":type(error).__name__,"reason":str(error)})
        save(output/"campaign_results.json",{"records":records,"manifest_sha256":digest(path),"prepare_only":prepare_only})
    journal.finalize()
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--manifest",type=Path)
    parser.add_argument("--record",action="append")
    parser.add_argument("--prepare-only",action="store_true")
    args=parser.parse_args()
    result=run(args.campaign_dir,args.output,manifest_path=args.manifest,record_ids=args.record,prepare_only=args.prepare_only)
    print(json.dumps({"jobs":len(result),"failed":sum(r["execution_status"]=="FAILED" for r in result)}))
    if any(r["execution_status"]=="FAILED" for r in result):raise SystemExit(1)


if __name__=="__main__":main()
