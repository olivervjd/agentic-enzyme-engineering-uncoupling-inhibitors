import unittest

from herbicide_desensitization_agent.app.backends.precomputed import PrecomputedComplexBackend
from herbicide_desensitization_agent.app.schemas.models import Ligand, Pose, Provenance, StructureModel, TargetProtein
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry


PROVENANCE = [Provenance("test://precomputed", "fixture", "synthetic")]


class PrecomputedBackendTests(unittest.TestCase):
    def setUp(self):
        self.entry = TargetRegistry.default().get("AT2G45300", "glyphosate")
        self.structures = [
            StructureModel(f"m{i}", self.entry.agi, "test", 0.8, list(self.entry.required_context), PROVENANCE)
            for i in (1, 2)
        ]
        self.poses = [Pose("p1", "m1", "glyphosate", 0.8, [10], PROVENANCE)]
        self.backend = PrecomputedComplexBackend(self.entry.agi, self.structures, {"glyphosate": self.poses})

    def test_replays_only_matching_pose_ensemble(self):
        ligand = Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE)
        self.assertEqual(self.backend.dock(self.structures, ligand), self.poses)

    def test_rejects_target_mismatch(self):
        target = TargetProtein("AT3G48560", "wrong", "MALW", provenance=PROVENANCE)
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.backend.predict_ensemble(target, self.entry)


if __name__ == "__main__":
    unittest.main()
