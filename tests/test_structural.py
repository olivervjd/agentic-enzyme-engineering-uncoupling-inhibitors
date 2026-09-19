import importlib.util
import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.app.backends.structural import TMAlignStructuralMatcher


AVAILABLE = all(importlib.util.find_spec(module) for module in ("Bio", "numpy", "tmtools"))


@unittest.skipUnless(AVAILABLE, "Install the structure extra for coordinate tests")
class StructuralTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.matcher = TMAlignStructuralMatcher()

    def pdb(self, name, *, translation=0, mutated=True, moved_site=0, omit=None, chain="A"):
        lines = []
        serial = 1
        for position in range(1, 9):
            if position == omit:
                continue
            residue = "VAL" if position == 2 and mutated else "ALA"
            for atom, delta in (("N", -0.7), ("CA", 0.0), ("C", 0.7), ("O", 1.1)):
                x = position * 2.1 + delta + translation
                y = (position % 3) * 1.2 + (moved_site if position == 2 else 0)
                z = (position % 2) * 1.7
                lines.append(f"ATOM  {serial:5d} {atom:^4s} {residue} {chain}{position:4d}    "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {atom[0]}\n")
                serial += 1
        path = self.root / name
        path.write_text("".join(lines) + "END\n")
        return path

    def test_rigid_transform_and_mutation_identity(self):
        ref = self.pdb("ref.pdb", mutated=False)
        mutant = self.pdb("mutant.pdb", translation=20)
        metrics = self.matcher.compare(mutant, ref, target_sequence="AAAAAAAA", mutation="A2V", active_site_residues=[2])
        self.assertAlmostEqual(metrics["query_tm_score"], 1.0, places=5)
        self.assertAlmostEqual(metrics["alignment_lddt"], 1.0, places=5)
        self.assertLess(metrics["active_site_rmsd_angstrom"], 1e-5)
        self.assertLess(metrics["global_ca_rmsd_angstrom"], 1e-5)

    def test_wild_type_cannot_masquerade_as_mutant(self):
        ref = self.pdb("ref.pdb", mutated=False)
        with self.assertRaisesRegex(ValueError, "sequence/numbering"):
            self.matcher.compare(ref, ref, target_sequence="AAAAAAAA", mutation="A2V")

    def test_site_motion_survives_global_fit(self):
        ref = self.pdb("ref.pdb", mutated=False)
        mutant = self.pdb("mutant.pdb", moved_site=4)
        result = self.matcher.compare(mutant, ref, target_sequence="AAAAAAAA", mutation="A2V", active_site_residues=[2])
        self.assertGreater(result["active_site_rmsd_angstrom"], 2)
        self.assertLess(result["alignment_lddt"], 1)

    def test_missing_site_is_not_zero_rmsd(self):
        ref = self.pdb("ref.pdb", mutated=False)
        mutant = self.pdb("mutant.pdb", omit=3)
        result = self.matcher.compare(mutant, ref, target_sequence="AAAAAAAA", mutation="A2V", active_site_residues=[3])
        self.assertIsNone(result["active_site_rmsd_angstrom"])
        self.assertEqual(result["alignment_coverage"], 7 / 8)

    def test_wrong_chain_and_offset_are_rejected(self):
        ref = self.pdb("ref.pdb", mutated=False)
        mutant = self.pdb("mutant.pdb")
        with self.assertRaisesRegex(ValueError, "Chain"):
            self.matcher.compare(mutant, ref, mutant_chain="B")
        with self.assertRaises(ValueError):
            self.matcher.compare(mutant, ref, mutant_offset=76, target_sequence="AAAAAAAA", mutation="A2V")

    def test_missing_mutation_is_rejected(self):
        ref = self.pdb("ref.pdb", mutated=False)
        mutant = self.pdb("mutant.pdb", omit=2)
        with self.assertRaisesRegex(ValueError, "absent"):
            self.matcher.compare(mutant, ref, target_sequence="AAAAAAAA", mutation="A2V")


if __name__ == "__main__":
    unittest.main()
