from __future__ import annotations

from pathlib import Path


class TMAlignStructuralMatcher:
    """TM-align plus residue-mapped local geometry for single-substitution models.

    Chain IDs and numbering offsets are explicit. RMSD uses one global CA fit;
    the active site is never refitted independently to hide pocket motion.
    """

    method = "TM-align (tmtools); Biopython Kabsch fit; mapped CA lDDT; active-site N/CA/C/O RMSD"

    def compare(
        self, mutant, reference, *, mutant_chain="A", reference_chain="A",
        mutant_offset=0, reference_offset=0, active_site_residues=(),
        target_sequence=None, mutation=None, sequence_start=1, sequence_end=None,
    ):
        import numpy as np
        from Bio.SVDSuperimposer import SVDSuperimposer
        from tmtools import tm_align

        query = self._residues(mutant, mutant_chain, mutant_offset)
        ref = self._residues(reference, reference_chain, reference_offset)
        if target_sequence is not None:
            sequence_end = sequence_end or len(target_sequence)
            if not 1 <= sequence_start <= sequence_end <= len(target_sequence):
                raise ValueError("Invalid modeled sequence domain")
            self._validate_sequence(ref, target_sequence, None, sequence_start, sequence_end)
            self._validate_sequence(query, target_sequence, mutation, sequence_start, sequence_end)
        common = sorted(query.keys() & ref.keys())
        if len(common) < 3:
            raise ValueError("At least three matched CA residues are required")
        qxyz = np.array([query[p]["CA"] for p in sorted(query)], dtype=float)
        rxyz = np.array([ref[p]["CA"] for p in sorted(ref)], dtype=float)
        aligned = tm_align(qxyz, rxyz, "".join(query[p]["aa"] for p in sorted(query)),
                           "".join(ref[p]["aa"] for p in sorted(ref)))
        qmatched = np.array([query[p]["CA"] for p in common], dtype=float)
        rmatched = np.array([ref[p]["CA"] for p in common], dtype=float)
        fit = SVDSuperimposer()
        fit.set(rmatched, qmatched)
        fit.run()
        rotation, translation = fit.get_rotran()

        # lDDT is per-residue distance agreement, independent of rigid-body fitting.
        rd = np.linalg.norm(rmatched[:, None] - rmatched[None, :], axis=-1)
        qd = np.linalg.norm(qmatched[:, None] - qmatched[None, :], axis=-1)
        neighbours = (rd < 15.0) & ~np.eye(len(common), dtype=bool)
        difference = np.abs(rd - qd)
        agreement = sum((difference < cutoff).astype(float) for cutoff in (0.5, 1.0, 2.0, 4.0)) / 4.0
        counts = neighbours.sum(axis=1)
        valid = counts > 0
        lddt = float(((agreement * neighbours).sum(axis=1)[valid] / counts[valid]).mean()) if valid.any() else None

        site = sorted(set(active_site_residues))
        site_rmsd = None
        if site and all(p in query and p in ref and all(a in query[p] and a in ref[p]
                       for a in ("N", "CA", "C", "O")) for p in site):
            qa = np.array([query[p][a] for p in site for a in ("N", "CA", "C", "O")])
            ra = np.array([ref[p][a] for p in site for a in ("N", "CA", "C", "O")])
            site_rmsd = float(np.sqrt(np.mean(np.sum((qa @ rotation + translation - ra) ** 2, axis=1))))
        denominators = [len(query), len(ref)]
        if target_sequence is not None:
            denominators.append(sequence_end - sequence_start + 1)
        return {
            "query_tm_score": float(aligned.tm_norm_chain1),
            "target_tm_score": float(aligned.tm_norm_chain2),
            "alignment_lddt": lddt,
            "alignment_coverage": min(len(common) / length for length in denominators),
            "global_ca_rmsd_angstrom": float(fit.get_rms()),
            "active_site_rmsd_angstrom": site_rmsd,
        }

    @staticmethod
    def _residues(path, chain, offset):
        import numpy as np
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.SeqUtils import seq1

        path = Path(path)
        parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
        model = parser.get_structure("model", str(path))[0]
        if chain not in model:
            raise ValueError(f"Chain {chain!r} not found in {path}")
        residues = {}
        for residue in model[chain]:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            if residue.id[2] != " ":
                raise ValueError("Insertion codes require an explicit residue map; integer offsets are insufficient")
            position = residue.id[1] + offset
            if position in residues:
                raise ValueError("Ambiguous residue numbering")
            atoms = {atom.name: np.asarray(atom.coord, dtype=float) for atom in residue}
            if any(not np.isfinite(x).all() for x in atoms.values()):
                raise ValueError("Non-finite structure coordinates")
            residues[position] = {"aa": seq1(residue.resname), **atoms}
        if len(residues) < 3:
            raise ValueError("Selected chain contains fewer than three protein CA atoms")
        return residues

    @staticmethod
    def _validate_sequence(residues, sequence, mutation, start, end):
        from ..validators.mutation_validator import validate_mutation

        expected = sequence
        if mutation and mutation != "WT":
            validate_mutation(mutation, sequence, set())
            position = int(mutation[1:-1])
            if position not in residues:
                raise ValueError("Mutated residue is absent from the selected structural chain")
            expected = sequence[:position - 1] + mutation[-1] + sequence[position:]
        for position, residue in residues.items():
            if not start <= position <= end or residue["aa"] != expected[position - 1]:
                raise ValueError(f"Structure sequence/numbering mismatch at position {position}")
