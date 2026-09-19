from __future__ import annotations


ADVERSARIAL_CASES = {
    "incorrect_ligand_pose": "reject-low-pose-confidence",
    "missing_cofactors": "reject-missing-required-context",
    "missing_membrane_or_complex": "reject-structure-context",
    "apo_structure_used_improperly": "reject-missing-native-ligand",
    "misleading_residue_numbering": "reject-sequence-structure-mapping",
    "essential_catalytic_residue": "reject-protected-residue",
    "shared_herbicide_native_contact": "reject-shared-contact",
    "conflicting_literature": "request-more-evidence",
    "low_confidence_structure": "request-better-structure",
    "tir1_conventional_inhibitor_assumption": "reject-mechanism-mismatch",
}


def evaluate_adversarial_signals(signals: dict[str, bool]) -> dict[str, str]:
    unknown = set(signals) - set(ADVERSARIAL_CASES)
    if unknown:
        raise ValueError(f"Unknown adversarial cases: {sorted(unknown)}")
    return {name: ADVERSARIAL_CASES[name] for name, detected in signals.items() if detected}
