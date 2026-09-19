import unittest

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
    def test_missing_direct_evidence_never_passes(self):
        record = FunctionRetentionAgent().assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[0]
        self.assertEqual(record.decision, "INSUFFICIENT_EVIDENCE")
        self.assertIsNone(record.ligand_kd_molar["glyphosate"])

    def test_direct_kd_and_structural_evidence_can_meet_screen(self):
        backend = PrecomputedFunctionRetentionBackend({"A2V": {
            "structural_metrics": {
                "query_tm_score": 0.95, "target_tm_score": 0.95, "alignment_lddt": 0.90,
                "alignment_coverage": 0.98, "active_site_rmsd_angstrom": 0.6,
            },
            "ligand_kd_molar": {"glyphosate": 2e-4, "phosphoenolpyruvate": 2e-6},
            "ligand_kd_fold_change_vs_wt": {"glyphosate": 20.0, "phosphoenolpyruvate": 1.5},
            "fold_ddg_kcal_mol": 0.8,
            "provenance": PROVENANCE,
        }})
        record = FunctionRetentionAgent(backend).assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[0]
        self.assertEqual(record.decision, "MEETS_COMPUTATIONAL_SCREEN")
        self.assertEqual(record.missing_evidence, [])

    def test_pic50_is_not_accepted_as_kd(self):
        backend = PrecomputedFunctionRetentionBackend({"A2V": {"affinity_pic50": 6.0}})
        record = FunctionRetentionAgent(backend).assess(
            TargetProtein("AT2G45300", "EPSPS", "MAV", provenance=PROVENANCE),
            [candidate()], [score()], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE),
            [Ligand("phosphoenolpyruvate", "native", "C", provenance=PROVENANCE)], [],
        )[0]
        self.assertEqual(record.decision, "INSUFFICIENT_EVIDENCE")


if __name__ == "__main__":
    unittest.main()
