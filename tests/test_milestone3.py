import unittest

from herbicide_desensitization_agent.app.agents.mutation_scoring import EvidenceAwareScoringAgent, pareto_rank
from herbicide_desensitization_agent.app.agents.pipeline_agents import ConstrainedMutationAgent
from herbicide_desensitization_agent.app.schemas.models import (
    InteractionFingerprint, MutationCandidate, Pose, Provenance, ScorePacket, TargetProtein,
)


PROVENANCE = [Provenance("test://m3", "fixture", "synthetic")]


class Milestone3Tests(unittest.TestCase):
    def test_generator_uses_only_fingerprint_guided_unprotected_positions(self):
        fingerprint = InteractionFingerprint("AT2G45300", [3], [4], [5], [1, 2], [6], PROVENANCE)
        target = TargetProtein("AT2G45300", "EPSPS", "MALWMRLL", provenance=PROVENANCE)
        candidates, rejected = ConstrainedMutationAgent().propose(target, {1, 2}, fingerprint)
        self.assertEqual({item.structure_residue for item in candidates}, {3})
        expanded, _ = ConstrainedMutationAgent(include_second_shell=True).propose(target, {1, 2}, fingerprint)
        self.assertEqual({item.structure_residue for item in expanded}, {3, 6})
        self.assertTrue(all(len(item.mutation) >= 3 for item in candidates))
        self.assertEqual(rejected, [])

    def test_scoring_preserves_components_evidence_and_uncertainty(self):
        fingerprint = InteractionFingerprint("AT2G45300", [3], [], [], [], [], PROVENANCE)
        candidate = MutationCandidate("L3I", "AT2G45300", 3, "test", "HERBICIDE_SELECTIVE_MUTABLE", PROVENANCE, 3)
        native = [Pose("n", "m", "native", 0.8, [7], PROVENANCE)]
        herbicide = [Pose("h", "m", "herbicide", 0.8, [3], PROVENANCE)]
        score = EvidenceAwareScoringAgent().score(candidate, native, herbicide, fingerprint)
        self.assertGreater(score.herbicide_escape_score, score.conservation_score)
        self.assertEqual(len(score.component_evidence), 9)
        self.assertTrue(any("binding free energies" in item for item in score.uncertainty))

    def test_pareto_ranking_does_not_collapse_components(self):
        candidate_a = MutationCandidate("A1V", "AT2G45300", 1, "a", "SECOND_SHELL_CANDIDATE", PROVENANCE)
        candidate_b = MutationCandidate("A2V", "AT2G45300", 2, "b", "SECOND_SHELL_CANDIDATE", PROVENANCE)
        candidate_c = MutationCandidate("A3V", "AT2G45300", 3, "c", "SECOND_SHELL_CANDIDATE", PROVENANCE)
        scores = [
            ScorePacket(*([0.9, 0.4] + [0.7] * 7), provenance=PROVENANCE),
            ScorePacket(*([0.4, 0.9] + [0.7] * 7), provenance=PROVENANCE),
            ScorePacket(*([0.3, 0.3] + [0.6] * 7), provenance=PROVENANCE),
        ]
        ranking = {item.mutation: item.front for item in pareto_rank([candidate_a, candidate_b, candidate_c], scores)}
        self.assertEqual(ranking["A1V"], 1)
        self.assertEqual(ranking["A2V"], 1)
        self.assertEqual(ranking["A3V"], 2)


if __name__ == "__main__":
    unittest.main()
