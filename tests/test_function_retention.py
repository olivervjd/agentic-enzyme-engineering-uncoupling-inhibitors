import unittest
from copy import deepcopy

from herbicide_desensitization_agent.app.agents.function_retention import FunctionRetentionAgent
from herbicide_desensitization_agent.app.backends.precomputed import PrecomputedFunctionRetentionBackend
from herbicide_desensitization_agent.app.schemas.models import (
    Ligand, MutationCandidate, Provenance, ScorePacket, TargetProtein,
)


PROVENANCE = [Provenance("test://retention", "fixture", "synthetic")]


def candidate():
    return MutationCandidate("A2V", "AT2G45300", 2, "test", "SECOND_SHELL_CANDIDATE", PROVENANCE, 2)


def score():
    return ScorePacket(*([0.8] * 9), provenance=PROVENANCE)


class FunctionRetentionTests(unittest.TestCase):
    def evidence(self):
        common = {
            "target_agi": "AT2G45300", "provenance": PROVENANCE,
            "affinity_protocol": {"glyphosate": "fixture-v1", "phosphoenolpyruvate": "fixture-v1"},
            "interval_description": "Synthetic fixture bounds, not scientific predictions",
            "structural_method": "synthetic structural comparison",
            "structural_metrics": {
                "query_tm_score": 0.99, "target_tm_score": 0.99, "alignment_lddt": 0.98,
                "alignment_coverage": 1.0, "active_site_rmsd_angstrom": 0.2, "global_ca_rmsd_angstrom": 0.3,
            },
            "fold_ddg_kcal_mol": 0.8,
        }
        wt, mutant = deepcopy(common), deepcopy(common)
        wt["ligand_kd_molar"] = {"glyphosate": 1e-5, "phosphoenolpyruvate": 1e-6}
        wt["ligand_kd_intervals_molar"] = {"glyphosate": [0.9e-5, 1.1e-5], "phosphoenolpyruvate": [0.9e-6, 1.1e-6]}
        mutant["ligand_kd_molar"] = {"glyphosate": 2e-4, "phosphoenolpyruvate": 1.5e-6}
        mutant["ligand_kd_intervals_molar"] = {"glyphosate": [1.8e-4, 2.2e-4], "phosphoenolpyruvate": [1.3e-6, 1.7e-6]}
        return {"WT": wt, "A2V": mutant}

    def assess(self, evidence):
        return FunctionRetentionAgent(PrecomputedFunctionRetentionBackend(evidence)).assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )

    def test_missing_direct_evidence_never_passes(self):
        record = FunctionRetentionAgent().assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[1]
        self.assertEqual(record.decision, "INSUFFICIENT_EVIDENCE")
        self.assertIsNone(record.ligand_kd_molar["glyphosate"])

    def test_direct_kd_and_structural_evidence_can_meet_screen(self):
        backend = PrecomputedFunctionRetentionBackend(self.evidence())
        record = FunctionRetentionAgent(backend).assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[1]
        self.assertEqual(record.decision, "MEETS_COMPUTATIONAL_SCREEN")
        self.assertEqual(record.missing_evidence, [])

    def test_pic50_is_not_accepted_as_kd(self):
        backend = PrecomputedFunctionRetentionBackend({"A2V": {"target_agi": "AT2G45300", "affinity_pic50": 6.0}})
        record = FunctionRetentionAgent(backend).assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[1]
        self.assertEqual(record.decision, "INSUFFICIENT_EVIDENCE")

    def test_wt_reference_and_ratios_are_computed(self):
        wt, mutant = self.assess(self.evidence())
        self.assertEqual(wt.mutation, "WT")
        self.assertEqual(wt.decision, "WILD_TYPE_REFERENCE")
        self.assertAlmostEqual(mutant.ligand_kd_fold_change_vs_wt["glyphosate"], 20)

    def test_native_stronger_binding_is_also_rejected(self):
        data = self.evidence()
        data["A2V"]["ligand_kd_molar"]["phosphoenolpyruvate"] = 1e-9
        data["A2V"]["ligand_kd_intervals_molar"]["phosphoenolpyruvate"] = [0.8e-9, 1.2e-9]
        self.assertIn("native_equivalence.phosphoenolpyruvate", self.assess(data)[1].failed_checks)

    def test_wide_uncertainty_is_not_a_pass(self):
        data = self.evidence()
        data["A2V"]["ligand_kd_intervals_molar"]["phosphoenolpyruvate"] = [0.1e-6, 10e-6]
        self.assertEqual(self.assess(data)[1].decision, "FAILS_COMPUTATIONAL_SCREEN")

    def test_missing_wt_and_mismatched_protocol_block_pass(self):
        for field in ("WT", "protocol"):
            data = self.evidence()
            if field == "WT":
                del data["WT"]
            else:
                data["WT"]["affinity_protocol"]["glyphosate"] = "different-model"
            self.assertEqual(self.assess(data)[1].decision, "INSUFFICIENT_EVIDENCE")

    def test_forged_ratios_are_rejected(self):
        data = self.evidence()
        data["A2V"]["ligand_kd_fold_change_vs_wt"] = {"glyphosate": 100}
        with self.assertRaisesRegex(ValueError, "disagrees"):
            self.assess(data)

    def test_nonphysical_values_are_rejected(self):
        for value in (float("nan"), float("inf"), -1, 0, "1e-6", True):
            with self.subTest(value=value):
                data = self.evidence()
                data["A2V"]["ligand_kd_molar"]["glyphosate"] = value
                with self.assertRaises(ValueError):
                    self.assess(data)

    def test_missing_intervals_never_pass(self):
        data = self.evidence()
        del data["A2V"]["ligand_kd_intervals_molar"]
        self.assertEqual(self.assess(data)[1].decision, "INSUFFICIENT_EVIDENCE")

    def test_target_mismatch_is_rejected(self):
        data = self.evidence()
        data["A2V"]["target_agi"] = "AT3G62980"
        with self.assertRaisesRegex(ValueError, "target_agi"):
            self.assess(data)

    def test_report_has_reference_units_and_legend(self):
        from herbicide_desensitization_agent.app.agents.function_retention import retention_csv, retention_markdown
        records = self.assess(self.evidence())
        markdown = retention_markdown(records, "glyphosate", ["phosphoenolpyruvate"])
        self.assertIn("| WT |", markdown)
        self.assertIn("Table legend:", markdown)
        self.assertIn("glyphosate Kd (M)", markdown)
        self.assertIn("phosphoenolpyruvate_ratio_lower", retention_csv(records, "glyphosate", ["phosphoenolpyruvate"]))


if __name__ == "__main__":
    unittest.main()
