from __future__ import annotations

from math import dist


def _pdb_atoms(pdb: str) -> list[tuple[int, tuple[float, float, float]]]:
    atoms = []
    for line in pdb.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        try:
            residue = int(line[22:26])
            coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atoms.append((residue, coordinates))
    return atoms


def _sdf_atoms(sdf: str) -> list[tuple[float, float, float]]:
    lines = sdf.splitlines()
    if len(lines) < 4:
        return []
    try:
        atom_count = int(lines[3][:3])
    except ValueError:
        return []
    atoms = []
    for line in lines[4 : 4 + atom_count]:
        try:
            atoms.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
        except ValueError:
            continue
    return atoms


def residue_contacts_from_pdb_and_sdf(pdb: str, sdf: str, cutoff_angstrom: float = 4.5) -> list[int]:
    if cutoff_angstrom <= 0:
        raise ValueError("contact cutoff must be positive")
    protein_atoms = _pdb_atoms(pdb)
    ligand_atoms = _sdf_atoms(sdf)
    return sorted(
        {
            residue
            for residue, protein_coordinates in protein_atoms
            if any(dist(protein_coordinates, ligand_coordinates) <= cutoff_angstrom for ligand_coordinates in ligand_atoms)
        }
    )

