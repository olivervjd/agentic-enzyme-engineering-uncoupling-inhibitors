from __future__ import annotations

from math import dist, isfinite


def _pdb_atoms(pdb: str, protein_chain: str | None = None) -> list[tuple[int, tuple[float, float, float]]]:
    atoms = []
    chains = {line[21] for line in pdb.splitlines() if line.startswith("ATOM  ") and len(line) >= 54}
    if protein_chain is None and len(chains) > 1:
        raise ValueError("Select protein_chain explicitly for a multichain structure")
    if protein_chain is not None and protein_chain not in chains:
        raise ValueError("Selected protein chain is absent")
    for line in pdb.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        if protein_chain is not None and line[21] != protein_chain:
            continue
        if line[16] not in {" ", "A"}:
            continue
        if line[26] != " ":
            raise ValueError("Insertion codes require an explicit residue map")
        element = line[76:78].strip() or line[12:16].strip().lstrip("0123456789")[:1]
        if element in {"H", "D"}:
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


def residue_contacts_from_pdb_and_sdf(
    pdb: str, sdf: str, cutoff_angstrom: float = 4.5,
    protein_chain: str | None = None, residue_offset: int = 0,
) -> list[int]:
    if not isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise ValueError("contact cutoff must be positive")
    protein_atoms = _pdb_atoms(pdb, protein_chain)
    ligand_atoms = _sdf_atoms(sdf)
    if not protein_atoms or not ligand_atoms:
        raise ValueError("Readable protein PDB and ligand V2000 SDF coordinates are required")
    return sorted(
        {
            residue + residue_offset
            for residue, protein_coordinates in protein_atoms
            if any(dist(protein_coordinates, ligand_coordinates) <= cutoff_angstrom for ligand_coordinates in ligand_atoms)
        }
    )
