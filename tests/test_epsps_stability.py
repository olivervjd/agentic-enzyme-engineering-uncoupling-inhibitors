"""Synthetic parser/protocol safeguards; actual DDGun outputs are run artifacts."""
import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.examples.epsps_stability import (
    read_a3m, map_mutations, compatibility_profile_module, parse_predictions,
    KNOWN_CONTROLS, HISTORICAL_HYPOTHESES,
)


class StabilityProtocolTests(unittest.TestCase):
    def sequence(self):
        sequence = list("A" * 444)
        for mutation in KNOWN_CONTROLS + HISTORICAL_HYPOTHESES:
            sequence[int(mutation[1:-1]) - 77] = mutation[0]
        return "".join(sequence)

    def test_a3m_insertions_removed_without_inventing_alignment(self):
        rows = read_a3m(">query\nACDE\n>homolog\nAcqC.D-\n")
        self.assertEqual(rows, [("query", "ACDE"), ("homolog", "ACD-")])
        with self.assertRaises(ValueError):
            read_a3m(">query\nACDE\n>wrong_width\nACD\n")
        with self.assertRaises(ValueError):
            read_a3m(">query\nACDE\n>query_again\nACDE\n")

    def test_all_mutations_match_exact_canonical_offset(self):
        rows = map_mutations(self.sequence())
        self.assertEqual(len(rows), 9)
        self.assertEqual(rows[0]["normalized_mutation"], "G101A")
        self.assertEqual(rows[1]["normalized_mutation"], "T102I")
        self.assertEqual(rows[2]["normalized_mutation"], "P106S")
        self.assertEqual({r["scope"] for r in rows[:3]}, {"known_control"})
        self.assertEqual({r["scope"] for r in rows[3:]}, {"historical_exploratory_hypothesis"})
        self.assertTrue(all(r["new_recommendation"] is False for r in rows))
        with self.assertRaises(ValueError):
            map_mutations(self.sequence(), ["A177V"])
        with self.assertRaises(ValueError):
            map_mutations(self.sequence(), ["G177A"], offset=75)

    def test_only_documented_compatibility_changes(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            (root / "tools").mkdir()
            text = ('try:\n        from Bio.SearchIO._legacy import NCBIStandalone\nexcept:\n        from Bio.Blast import NCBIStandalone\n'
                    'handle = open(alifile, "rU")\nscore = 0.18 * feature\n')
            original = root / "tools" / "ali2prof.py"
            original.write_text(text)
            result = compatibility_profile_module(root, root / "compatibility")
            self.assertEqual(original.read_text(), text)
            patched = Path(result["module"]).read_text()
            self.assertIn('open(alifile, "r")', patched)
            self.assertIn("score = 0.18 * feature", patched)
            self.assertNotIn("from Bio", patched)
            self.assertFalse(result["algorithm_changes"])
            self.assertNotEqual(result["original_sha256"], result["patched_sha256"])

    def test_native_output_sign_reversal_and_missing_uncertainty(self):
        mapping = map_mutations(self.sequence(), ["G177A", "T178I"])
        text = "#PDBFILE\tCHAIN\tVARIANT\tS_DDG[3D]\tT_DDG[3D]\tSTABILITY[3D]\nfile.pdb\tA\tG101A\t0.5\t0.5\tIncrease\nfile.pdb\tA\tT102I\t-0.7\t-0.7\tDecrease\n"
        rows = parse_predictions(text, mapping)
        self.assertEqual(rows[0]["folding_ddg_kcal_mol"], -.5)
        self.assertEqual(rows[1]["folding_ddg_kcal_mol"], .7)
        self.assertIsNone(rows[0]["uncertainty_interval"])
        self.assertEqual(rows[0]["calibration_status"], "UNCALIBRATED_CUSTOM_MSA_PROTOCOL")
        with self.assertRaises(ValueError):
            parse_predictions(text.replace("Increase", "Decrease"), mapping)
        with self.assertRaises(ValueError):
            parse_predictions(text.split("file.pdb\tA\tT102I")[0], mapping)


if __name__ == "__main__":
    unittest.main()
