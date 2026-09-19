import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.app.agents.learning_agent import LearningAgent
from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.orchestrator.learning_workflow import LearningWorkflow
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import AssayResult, Provenance
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


PROVENANCE = [Provenance("test://assay", "synthetic-fixture", "synthetic")]


class Milestone5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = TargetRegistry.default()
        cls.entry = cls.registry.get("AT2G45300", "glyphosate")
        backend = MockScientificBackend()
        cls.workflow = WorkflowOrchestrator(cls.registry, backend, backend, backend, backend, backend)
        cls.request = make_request(cls.entry)
        cls.result = cls.workflow.run(cls.request)
        cls.mutation = cls.result.packets[0].candidate.mutation

    def assay(self, **measurement_changes):
        measurements = {
            "herbicide_response_fold_change": 3.0,
            "protein_expression_fraction": 0.9,
            "native_activity_fraction": 0.8,
            "thermal_stability_delta_c": -1.0,
        }
        measurements.update(measurement_changes)
        units = {name: "fraction" for name in measurements}
        return AssayResult(
            "assay-1", self.entry.agi, self.entry.herbicide, self.mutation,
            measurements, units, PROVENANCE,
        )

    def test_recalibration_is_component_wise_and_keeps_warnings(self):
        report = LearningAgent(self.registry).recalibrate(
            self.entry.agi, self.entry.herbicide, self.result.packets,
            [self.assay()], self.request.target.sequence,
        )
        self.assertEqual(report.assay_count, 1)
        self.assertIn("herbicide_escape_score", report.component_bias)
        self.assertIn(self.mutation, report.calibrated_components)
        self.assertTrue(report.warnings)

    def test_next_round_excludes_measured_mutation_and_requires_review(self):
        agent = LearningAgent(self.registry)
        assay = self.assay()
        report = agent.recalibrate(
            self.entry.agi, self.entry.herbicide, self.result.packets, [assay], self.request.target.sequence
        )
        selected = agent.select_next_round(self.result.packets, [assay], report)
        self.assertNotIn(self.mutation, {item.mutation for item in selected})
        self.assertTrue(all(item.status == "NEEDS_REVIEW" for item in selected))

    def test_missing_mechanism_specific_measurement_is_rejected(self):
        assay = self.assay()
        measurements = dict(assay.measurements)
        measurements.pop("native_activity_fraction")
        incomplete = AssayResult(
            assay.assay_id, assay.agi, assay.herbicide, assay.mutation,
            measurements, {name: "fraction" for name in measurements}, assay.provenance,
        )
        with self.assertRaisesRegex(ValueError, "missing measurements"):
            LearningAgent(self.registry).recalibrate(
                self.entry.agi, self.entry.herbicide, self.result.packets,
                [incomplete], self.request.target.sequence,
            )

    def test_learning_workflow_writes_required_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            learning = LearningWorkflow(LearningAgent(self.registry), ArtifactStore(directory))
            learning.run(self.request, self.result, [self.assay()])
            run_dir = Path(directory) / "AT2G45300-glyphosate"
            self.assertTrue((run_dir / "assay_results.json").is_file())
            self.assertTrue((run_dir / "model_recalibration_report.json").is_file())
            self.assertTrue((run_dir / "next_round_candidates.json").is_file())

    def test_recalibration_requires_measurements(self):
        with self.assertRaisesRegex(ValueError, "At least one assay"):
            LearningAgent(self.registry).recalibrate(
                self.entry.agi, self.entry.herbicide, self.result.packets,
                [], self.request.target.sequence,
            )


if __name__ == "__main__":
    unittest.main()
