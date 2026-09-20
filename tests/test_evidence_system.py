"""Synthetic numerical fixtures verify software gates, never scientific performance."""
import copy
import json
import unittest
from unittest.mock import Mock

from herbicide_desensitization_agent.app.evidence_system import (
    assess_binding, affinity_difference, calibrate_thresholds, calibration_errors,
    classify_residues, compare_contexts, evaluate_candidate, extract_pose_contacts,
    summarize_contacts, validate_stage,
)
from herbicide_desensitization_agent.app.evidence_system.schema import STRUCTURAL_METRICS, digest

PROVENANCE = {"source_type": "computational", "source": "synthetic tool fixture", "artifact": "fixture://test-only"}


def context(mutant=False, method="method-a", seed=1):
    return {"protein_sequence": "ACD" if not mutant else "AAD", "substitutions": ["C2A"] if mutant else [],
            "isoform": "test", "modeled_residue_range": [1, 3], "residue_mapping": {"A:1": 1, "A:2": 2, "A:3": 3},
            "oligomeric_state": "monomer", "cellular_context": "soluble", "substrates": ["native"],
            "cofactors": [], "ligands": [{"id": "test", "canonical_identifier": "test-only", "stereochemistry": "achiral",
                                           "protonation": "specified", "tautomer": "specified", "formal_charge": 0}],
            "pH": 7.0, "ionic_assumptions": "test", "metal_ions": [], "catalytic_waters": [], "partners": [],
            "structure_model": method, "structure_model_version": "1", "affinity_model": method,
            "affinity_model_version": "1", "msa_source": "test-only", "msa_depth": 10, "seed": seed, "sampling": {"n": 3}}


def affinity(estimate, half=.01, metric="pIC50", unit="dimensionless"):
    return {"metric": metric, "unit": unit, "estimate": estimate, "interval": [estimate-half, estimate+half],
            "interval_method": "synthetic software fixture", "interval_level": .95, "provenance": PROVENANCE.copy()}


def calibration():
    metrics = {name: {"direction": "min" if name in ("tm_score_forward", "tm_score_reverse", "ca_lddt", "coverage") else "max",
                       "limit": .8 if name in ("tm_score_forward", "tm_score_reverse", "ca_lddt", "coverage") else 1,
                       "margin": .1, "unit": "dimensionless"} for name in STRUCTURAL_METRICS}
    spec = {"metrics": metrics, "assembly_required": False, "required_evaluations": {"herb": "herbicide", "native": "native_substrate"},
            "binding": {"minimum_replicates": 2, "minimum_methods": 2, "pose_consistency_min": .8, "contact_reproducibility_min": .8,
                        "evaluations": {"herb": {"metric": "pIC50", "unit": "dimensionless", "weakening_direction": "decrease",
                                                 "weakening_threshold": .5, "wt_variability_bound": .2, "binding_support_interval": [4, 10]},
                                        "native": {"metric": "pIC50", "unit": "dimensionless", "equivalence_interval": [-.2, .2],
                                                   "wt_variability_bound": .2, "binding_support_interval": [4, 10]}}}}
    values = {name: .95 if rule["direction"] == "min" else .05 for name, rule in metrics.items()}
    wt = [{"method": "a", "seed": i, "metrics": values, "provenance": PROVENANCE.copy(),
           "binding_measurements": {name: {"value": 5 + i*.01, "metric": "pIC50", "unit": "dimensionless"} for name in ("herb", "native")}} for i in (1,2,3)]
    controls = [{"id": kind, "kind": kind, "expected": "pass", "metrics": values, "provenance": PROVENANCE.copy()}
                for kind in ("resistant", "catalytic_loss", "positive_binding", "negative_binding")]
    for control in controls:
        kind = control["kind"]
        if kind in ("resistant", "catalytic_loss"):
            control["ligand_changes"] = {ligand: {"metric": "pIC50", "unit": "dimensionless",
                                                "value": -1 if (kind == "resistant" and ligand == "herb") or (kind == "catalytic_loss" and ligand == "native") else 0}
                                         for ligand in ("herb", "native")}
        else:
            control["ligand_affinities"] = {ligand: {"metric": "pIC50", "unit": "dimensionless", "value": 5 if kind == "positive_binding" else 1}
                                           for ligand in ("herb", "native")}
    return calibrate_thresholds(spec, wt, controls, registered_at="2026-01-01T00:00:00Z", candidates_started_at="2026-01-02T00:00:00Z")


def assessment(evaluation="herb"):
    methods = []
    for family in ("a", "b"):
        methods.append({"independent_family": family, "context": context(True, family), "wt_context": context(False, family),
                        "conclusion": "supports_weakening" if evaluation == "herb" else "supports_binding",
                        "replicates": [{"seed": i, "artifact": f"fixture://{family}/{i}"} for i in (1,2,3)],
                        "affinity": affinity(4 if evaluation == "herb" else 5), "wt_affinity": affinity(5),
                        "pose_consistency": .95, "contact_reproducibility": .95,
                        "control_comparison": {"passed": True, "positive_control_id": "positive_binding", "negative_control_id": "negative_binding", "provenance": PROVENANCE.copy()}})
    return {"variant": "C2A", "ligand": evaluation, "evaluation_id": evaluation, "methods": methods}


def candidate():
    value = {"variant": "C2A", "required_evaluations": ["herb", "native"], "bindings": [],
             "structural_metrics": {name: {"value": .95 if name in ("tm_score_forward", "tm_score_reverse", "ca_lddt", "coverage") else .05,
                                           "unit": "dimensionless", "provenance": PROVENANCE.copy()} for name in STRUCTURAL_METRICS}}
    for ligand, role in (("herb", "herbicide"), ("native", "native_substrate")):
        a = assessment(ligand)
        value["bindings"].append({"evaluation_id": ligand, "ligand": ligand, "role": role,
                                  "wt_context": context(False, "a"), "mutant_context": context(True, "a"),
                                  "wt_affinity": affinity(5), "mutant_affinity": a["methods"][0]["affinity"], "assessment": a})
    for name in ("biological_context_complete", "contacts_reproducible", "protected_interactions_retained"):
        value[name] = {"passed": True, "provenance": PROVENANCE.copy(), "reason": "Synthetic software fixture"}
    return value


class EvidenceSystemTests(unittest.TestCase):
    def test_complete_context_pair_and_cross_method(self):
        self.assertTrue(compare_contexts(context(), context(True), mutant=True)["equivalent"])
        self.assertTrue(compare_contexts(context(method="a"), context(method="b"), cross_method=True)["equivalent"])
        self.assertFalse(compare_contexts(context(), context(seed=2))["equivalent"])
        other = context(True)
        other["substitutions"] = ["C3A"]
        self.assertFalse(compare_contexts(context(), other, mutant=True)["equivalent"])

    def test_missing_unknown_and_chemical_mismatch(self):
        for key in ("pH", "msa_depth", "residue_mapping", "seed", "sampling"):
            other = context()
            other.pop(key)
            self.assertFalse(compare_contexts(other, other)["equivalent"])
        other = context()
        other["ligands"][0]["formal_charge"] = -1
        self.assertFalse(compare_contexts(context(), other)["equivalent"])

    def test_calibration_reproduces_and_tampering_rejected(self):
        cal = calibration()
        self.assertEqual(cal["status"], "VALIDATED", cal["limitations"])
        self.assertEqual(calibration_errors(cal), [])
        cal["thresholds"]["folding_ddg"]["limit"] = 100
        cal["digest"] = digest({k:v for k,v in cal.items() if k != "digest"})
        self.assertTrue(calibration_errors(cal))

    def test_no_posthoc_calibration_or_empty_controls(self):
        cal = calibration()
        with self.assertRaises(ValueError):
            calibrate_thresholds(cal["specification"], cal["wt_replicates"], cal["controls"], registered_at="2026-02-01T00:00:00Z", candidates_started_at="2026-01-01T00:00:00Z")
        spec = copy.deepcopy(cal["specification"])
        spec["required_control_kinds"] = []
        no_controls = calibrate_thresholds(spec, cal["wt_replicates"], [], registered_at=cal["registered_at"], candidates_started_at=cal["candidates_started_at"])
        self.assertNotEqual(no_controls["status"], "VALIDATED")

    def test_duplicate_wt_replicates_and_incorrect_direction_rejected(self):
        cal = calibration()
        cal["wt_replicates"][1]["seed"] = 1
        cal["specification"]["binding"]["evaluations"]["herb"]["weakening_direction"] = "increase"
        output = calibrate_thresholds(cal["specification"], cal["wt_replicates"], cal["input_controls"], registered_at=cal["registered_at"], candidates_started_at=cal["candidates_started_at"])
        self.assertNotEqual(output["status"], "VALIDATED")

    def test_seven_evidence_states_no_pose_binding_claim(self):
        cal = calibration()
        self.assertEqual(assess_binding({}, cal)["evidence_state"], "Not evaluated")
        self.assertEqual(assess_binding({"evaluated": True}, cal)["evidence_state"], "Insufficient evidence")
        record = assessment()
        self.assertEqual(assess_binding(record, cal)["evidence_state"], "Strong computational support for weakened binding")
        self.assertEqual(assess_binding(assessment("native"), cal)["evidence_state"], "Strong computational support for binding")
        record["methods"] = record["methods"][:1]
        record["methods"][0]["replicates"] = record["methods"][0]["replicates"][:1]
        self.assertEqual(assess_binding(record, cal)["evidence_state"], "Insufficient evidence")

    def test_disagreement_only_in_equivalent_context(self):
        record = assessment()
        record["methods"][1]["conclusion"] = "supports_binding"
        self.assertEqual(assess_binding(record, calibration())["evidence_state"], "Conflicting computational evidence")
        record["methods"][1]["context"]["pH"] = 9
        output = assess_binding(record, calibration())
        self.assertEqual(output["evidence_state"], "Insufficient evidence")
        self.assertEqual(output["method_agreement"], "non_equivalent")

    def test_weakening_cannot_be_asserted_without_wt_comparison(self):
        record = assessment()
        for method in record["methods"]:
            method.pop("wt_affinity")
        self.assertEqual(assess_binding(record, calibration())["evidence_state"], "Insufficient evidence")

    def test_affinity_quantities_never_converted(self):
        self.assertFalse(affinity_difference(affinity(5), affinity(5, metric="Kd", unit="nM"))["valid"])
        self.assertFalse(affinity_difference(affinity(5), affinity(5, metric="docking_confidence"))["valid"])
        invalid = affinity(5)
        invalid["interval"] = [float("nan"), 6]
        self.assertFalse(affinity_difference(affinity(5), invalid)["valid"])
        interval = affinity_difference(affinity(5, .2), affinity(4, .3))["interval"]
        self.assertAlmostEqual(interval[0], -1.5)
        self.assertAlmostEqual(interval[1], -.5)

    def test_complete_screen_and_json_contract(self):
        result = evaluate_candidate(candidate(), calibration())
        self.assertEqual(result["decision"], "MEETS_COMPUTATIONAL_SCREEN", result["limitations"])
        json.dumps(result, allow_nan=False)

    def test_absent_evidence_is_not_pass(self):
        self.assertEqual(evaluate_candidate({"variant": "legacy"})["decision"], "INSUFFICIENT_EVIDENCE")
        cal, case = calibration(), candidate()
        case["structural_metrics"].pop("folding_ddg")
        self.assertEqual(evaluate_candidate(case, cal)["decision"], "INSUFFICIENT_EVIDENCE")

    def test_missing_native_cannot_be_removed_from_required_registry(self):
        case = candidate()
        case["required_evaluations"] = ["herb"]
        case["bindings"] = case["bindings"][:1]
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "INSUFFICIENT_EVIDENCE")

    def test_affinity_interval_entirely_inside_equivalence(self):
        case = candidate()
        case["bindings"][1]["mutant_affinity"]["interval"] = [4.5, 5.5]
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "FAILS_COMPUTATIONAL_SCREEN")

    def test_other_ligand_evidence_cannot_be_reused(self):
        case = candidate()
        case["bindings"][0]["assessment"]["ligand"] = "unrelated"
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "INSUFFICIENT_EVIDENCE")

    def test_unrelated_wt_affinity_cannot_drive_comparison(self):
        case = candidate()
        case["bindings"][0]["wt_affinity"] = affinity(9)
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "INSUFFICIENT_EVIDENCE")

    def test_control_label_without_ligand_measurements_does_not_calibrate(self):
        cal = calibration()
        for control in cal["input_controls"]:
            control.pop("ligand_changes", None)
            control.pop("ligand_affinities", None)
        output = calibrate_thresholds(cal["specification"], cal["wt_replicates"], cal["input_controls"],
                                      registered_at=cal["registered_at"], candidates_started_at=cal["candidates_started_at"])
        self.assertNotEqual(output["status"], "VALIDATED")

    def test_negative_pocket_volume_change_is_not_automatically_safe(self):
        case = candidate()
        case["structural_metrics"]["pocket_volume_change"]["value"] = -100
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "FAILS_COMPUTATIONAL_SCREEN")

    def test_llm_measurement_cannot_pass(self):
        case = candidate()
        case["structural_metrics"]["folding_ddg"]["provenance"]["source_type"] = "language_model"
        self.assertEqual(evaluate_candidate(case, calibration())["decision"], "INSUFFICIENT_EVIDENCE")

    def test_experimental_binding_is_not_kinetic_or_text_claim(self):
        record = {"experimental": {"provenance": {**PROVENANCE, "source_type": "experimental"}, "context_matched": True,
                                   "outcome": "binding", "assay": "test", "direct_binding_assay": True, "measurement": "binds"}}
        self.assertNotEqual(assess_binding(record)["evidence_state"], "Experimentally confirmed binding")
        record["experimental"]["measurement"] = affinity(5)
        self.assertNotEqual(assess_binding(record)["evidence_state"], "Experimentally confirmed binding")
        measurement = affinity(5, metric="Kd", unit="nM")
        measurement["provenance"]["source_type"] = "experimental"
        record["experimental"]["measurement"] = measurement
        self.assertEqual(assess_binding(record)["evidence_state"], "Experimentally confirmed binding")

    def test_experimental_learning_does_not_autoapprove(self):
        case = candidate()
        case["experimental_validation"] = {"outcome": "validated", "reviewed": False, "scope": "herbicide_interference_and_native_function",
                                           "provenance": {**PROVENANCE, "source_type": "experimental"}, "measurement_artifacts": ["assay.csv"]}
        self.assertEqual(evaluate_candidate(case)["decision"], "INSUFFICIENT_EVIDENCE")
        case["experimental_validation"]["reviewed"] = True
        self.assertEqual(evaluate_candidate(case)["decision"], "EXPERIMENTALLY_VALIDATED")
        case["experimental_validation"]["outcome"] = "rejected"
        self.assertEqual(evaluate_candidate(case)["decision"], "EXPERIMENTALLY_REJECTED")

    def test_stage_validation(self):
        self.assertFalse(validate_stage("affinity", {"affinity": {"estimate": 1}})["valid"])
        self.assertTrue(validate_stage("affinity", {"affinity": affinity(5)})["valid"])


class ContactEvidenceTests(unittest.TestCase):
    def rows(self):
        protein = [{"name": "NZ", "element": "N", "xyz": [0,0,0], "chain": "A", "structure_residue": "12A", "amino_acid": "LYS"}]
        ligand = [{"name": "O1", "element": "O", "xyz": [3,0,0]}]
        return extract_pose_contacts(protein, ligand, variant="WT", ligand="herb", pose_id="p1", replicate_id="r1",
                                     model_id="m1", method="a", canonical_mapping={"A:12A": 15}, provenance=PROVENANCE)

    def test_coordinate_contacts_and_untyped_chemistry(self):
        row = self.rows()[0]
        self.assertEqual(row["canonical_residue"], 15)
        self.assertEqual(row["min_distance"], 3)
        self.assertEqual(row["interaction_types"], ["general_proximity"])

    def test_frequency_denominator_includes_absent_contacts(self):
        registry = [{"variant": "WT", "ligand": "herb", "pose_id": f"p{i}", "replicate_id": f"r{i}", "model_id": f"m{i}", "method": "a"} for i in (1,2)]
        summary = summarize_contacts(self.rows(), registry)[0]
        self.assertEqual(summary["contact_frequency"], .5)
        self.assertEqual(summary["independent_model_frequency"], .5)
        self.assertLess(summary["contact_interval"][0], .5)
        self.assertGreater(summary["contact_interval"][1], .5)

    def test_unregistered_contact_rejected(self):
        with self.assertRaises(ValueError):
            summarize_contacts(self.rows(), [])

    def test_differential_contacts_uncertainty_and_protection(self):
        thresholds = {"reproducibility": .7, "native_max": .2, "selectivity_min": .3, "method_tolerance": .2,
                      "conservation_max": .8, "minimum_models": 3, "minimum_methods": 2, "second_shell_max_angstrom": 6}
        strong = {"frequency": .95, "interval": [.8, 1], "method_frequencies": {"a": .95, "b": .95}, "n": 20}
        weak = {"frequency": .05, "interval": [0, .2], "method_frequencies": {"a": .05, "b": .05}, "n": 20}
        row = {"canonical_residue": 10, "protected": False, "catalytic": False, "conservation": .2, "cross_method_binding_site_support": True,
               "herbicide": strong, "native": {"native": weak}}
        def classify(value):
            return classify_residues([value], native_ligands=["native"], thresholds=thresholds)[0]
        self.assertEqual(classify(row)["classification"], "HERBICIDE_SELECTIVE_CONTACT")
        self.assertTrue(classify(row)["nomination_eligible"])
        self.assertEqual(classify({**row, "native": {}})["classification"], "INSUFFICIENT_EVIDENCE")
        self.assertEqual(classify({**row, "native": {"native": strong}})["classification"], "SHARED_HERBICIDE_NATIVE_CONTACT")
        self.assertEqual(classify({**row, "protected": True})["classification"], "CATALYTIC_OR_PROTECTED")
        self.assertFalse(classify({**row, "conservation": .99})["nomination_eligible"])
        self.assertEqual(classify({**row, "herbicide": weak, "sequence_adjacent": True})["classification"], "INSUFFICIENT_EVIDENCE")
        self.assertEqual(classify({**row, "herbicide": weak, "distance_to_herbicide_contact_angstrom": 4, "spatial_evidence": PROVENANCE})["classification"], "SECOND_SHELL_CANDIDATE")
        invalid = {**strong, "interval": [0, .2]}
        self.assertEqual(classify({**row, "herbicide": invalid})["classification"], "INSUFFICIENT_EVIDENCE")


class ProductionGuardTests(unittest.TestCase):
    def test_old_quality_pass_cannot_generate_live_mutations(self):
        from herbicide_desensitization_agent.app.agents.quality import WorkflowReassessmentRequired
        from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
        from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
        from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
        from herbicide_desensitization_agent.examples.run_all import make_request
        registry = TargetRegistry.default()
        backend, quality, evidence = MockScientificBackend(), Mock(), Mock()
        quality.assess.return_value = quality.assess_pocket.return_value = {"status": "PASSED", "reasons": []}
        evidence.synthesize.return_value = {"facts": []}
        workflow = WorkflowOrchestrator(registry, *([backend]*5), require_live_gates=True, quality_agent=quality, evidence_agent=evidence)
        workflow.mutation_agent = Mock()
        with self.assertRaises(WorkflowReassessmentRequired) as raised:
            workflow.run(make_request(registry.entries[0]))
        self.assertIn("calibration missing", str(raised.exception))
        workflow.mutation_agent.propose.assert_not_called()
        self.assertEqual(next(s["status"] for s in workflow.stage_manifest if s["stage"] == "calibrated_evidence_baseline"), "BLOCKED")


if __name__ == "__main__":
    unittest.main()
