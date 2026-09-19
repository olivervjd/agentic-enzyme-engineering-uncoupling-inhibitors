import unittest

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.evals.adversarial_cases import ADVERSARIAL_CASES, evaluate_adversarial_signals
from herbicide_desensitization_agent.app.evals.benchmark_runner import BenchmarkCase, FixedTargetBenchmarkRunner
from herbicide_desensitization_agent.app.evals.judges import RUBRIC, IndependentDeterministicJudge, RosalindDomainJudge
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.examples.run_all import make_request


class Milestone4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = TargetRegistry.default()
        cls.backend = MockScientificBackend()
        cls.workflow = WorkflowOrchestrator(
            cls.registry, cls.backend, cls.backend, cls.backend, cls.backend, cls.backend
        )
        cls.request = make_request(cls.registry.entries[0])
        cls.result = cls.workflow.run(cls.request)

    def test_two_independent_judgments_exist_per_packet(self):
        self.assertEqual(len(self.result.judge_results), 2 * len(self.result.packets))
        judge_ids = {item.judge_id for item in self.result.judge_results}
        self.assertEqual(judge_ids, {"mock-rosalind-domain-judge", "independent-deterministic-judge"})
        self.assertTrue(all(set(item.category_scores) == set(RUBRIC) for item in self.result.judge_results))

    def test_review_packet_has_required_sections(self):
        packet = self.result.packets[0]
        self.assertTrue(packet.mechanistic_hypothesis)
        self.assertTrue(packet.herbicide_interactions_disrupted)
        self.assertTrue(packet.native_function_interactions_preserved)
        self.assertEqual(set(packet.risk_summary), {"fold", "conservation", "cofactor_or_complex"})
        self.assertEqual(packet.status, "NEEDS_REVIEW")

    def test_fixed_target_benchmark_reports_all_metrics(self):
        case = BenchmarkCase(
            "synthetic-test", "AT2G45300", "glyphosate", self.request.target.sequence,
            self.request.protected_residues, {self.result.packets[0].candidate.structure_residue}, "synthetic",
        )
        report = FixedTargetBenchmarkRunner(self.registry).run(case, self.result.packets)
        self.assertEqual(report["protected_residue_violation_rate"], 0.0)
        self.assertEqual(report["invalid_mutation_rate"], 0.0)
        self.assertIn("top_5_known_position_recall", report)

    def test_benchmark_rejects_out_of_scope_pairing(self):
        case = BenchmarkCase("bad", "AT2G45300", "atrazine", "MALW", set(), {1}, "synthetic")
        with self.assertRaises(ValueError):
            FixedTargetBenchmarkRunner(self.registry).run(case, [])

    def test_adversarial_catalog_covers_required_cases(self):
        findings = evaluate_adversarial_signals({name: True for name in ADVERSARIAL_CASES})
        self.assertEqual(len(findings), 10)
        self.assertEqual(findings["tir1_conventional_inhibitor_assumption"], "reject-mechanism-mismatch")


if __name__ == "__main__":
    unittest.main()
