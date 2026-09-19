import unittest
from dataclasses import replace

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.evals.deterministic_checks import validate_evaluation_packet
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.validators.mutation_validator import validate_mutation
from herbicide_desensitization_agent.examples.run_all import make_request


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = TargetRegistry.default()
        backend = MockScientificBackend()
        cls.workflow = WorkflowOrchestrator(cls.registry, backend, backend, backend, backend, backend)

    def test_registry_contains_only_six_fixed_pairings(self):
        self.assertEqual(len(self.registry.entries), 6)
        self.assertEqual(
            {entry.agi for entry in self.registry.entries},
            {"AT2G45300", "AT1G66200", "AT5G35630", "ATCG00020", "AT3G48560", "AT3G62980"},
        )

    def test_every_pairing_runs_and_requires_review(self):
        for entry in self.registry.entries:
            with self.subTest(agi=entry.agi, herbicide=entry.herbicide):
                workflow_request = make_request(entry)
                result = self.workflow.run(workflow_request)
                self.assertEqual(result.packets[0].status, "NEEDS_MORE_COMPUTATION")
                self.assertEqual(result.function_retention[0].mutation, "WT")
                self.assertEqual(result.packets[0].provenance[0].evidence_type, "synthetic")
                self.assertEqual(
                    validate_evaluation_packet(
                        result.packets[0], workflow_request.target.sequence, workflow_request.protected_residues
                    ),
                    [],
                )

    def test_unknown_pairing_is_rejected(self):
        request = make_request(self.registry.entries[0])
        changed = request.__class__(
            target=request.target,
            herbicide=request.herbicide.__class__(
                "atrazine", "herbicide", "C", provenance=request.herbicide.provenance
            ),
            native_ligands=request.native_ligands,
            functional_context=request.functional_context,
            protected_residues=request.protected_residues,
        )
        with self.assertRaisesRegex(ValueError, "Unsupported AGI"):
            self.workflow.run(changed)

    def test_protected_mutation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Protected residue"):
            validate_mutation("M1L", "MALW", {1})

    def test_missing_native_pose_blocks_mutation_selection(self):
        class MissingNative(MockScientificBackend):
            def dock(self, structures, ligand):
                return [] if ligand.role == "native" else super().dock(structures, ligand)
        backend = MissingNative()
        workflow = WorkflowOrchestrator(self.registry, backend, backend, backend, backend, backend)
        with self.assertRaisesRegex(ValueError, "pose evidence"):
            workflow.run(make_request(self.registry.entries[0]))

    def test_glutamine_synthetase_requires_all_substrates(self):
        for agi in ("AT1G66200", "AT5G35630"):
            entry = self.registry.get(agi, "glufosinate")
            self.assertEqual(set(entry.native_ligands), {"glutamate", "ATP", "ammonium"})
            request = make_request(entry)
            request = replace(request, native_ligands=request.native_ligands[:1])
            with self.assertRaisesRegex(ValueError, "Missing native ligands"):
                self.workflow.run(request)

    def test_empty_ensemble_is_rejected(self):
        class EmptyStructures(MockScientificBackend):
            def predict_ensemble(self, target, entry):
                return []
        backend = EmptyStructures()
        workflow = WorkflowOrchestrator(self.registry, backend, backend, backend, backend, backend)
        with self.assertRaisesRegex(ValueError, "nonempty structure ensemble"):
            workflow.run(make_request(self.registry.entries[0]))


if __name__ == "__main__":
    unittest.main()
