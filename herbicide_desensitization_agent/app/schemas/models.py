from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


@dataclass(frozen=True)
class Provenance:
    source: str
    method: str
    evidence_type: str
    notes: str = ""


@dataclass(frozen=True)
class TargetRegistryEntry:
    herbicide: str
    native_ligands: list[str]
    protein_name: str
    symbol: str
    agi: str
    cellular_compartment: str
    mechanism_class: str
    required_context: list[str]
    modeling_workflow: str
    validation_rules: list[str]
    provenance: list[Provenance]


@dataclass(frozen=True)
class TargetProtein:
    agi: str
    name: str
    sequence: str
    residue_numbering_start: int = 1
    provenance: list[Provenance] = field(default_factory=list)


@dataclass(frozen=True)
class Ligand:
    name: str
    role: Literal["native", "herbicide", "cofactor"]
    structure: str
    structure_format: str = "SMILES"
    provenance: list[Provenance] = field(default_factory=list)


@dataclass(frozen=True)
class FunctionalPartner:
    name: str
    kind: str
    provenance: list[Provenance] = field(default_factory=list)


@dataclass(frozen=True)
class StructureModel:
    model_id: str
    target_agi: str
    method: str
    confidence: float
    context: list[str]
    provenance: list[Provenance]
    artifact_path: str | None = None
    format: str | None = None
    scores: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Pose:
    pose_id: str
    model_id: str
    ligand_name: str
    confidence: float
    contacts: list[int]
    provenance: list[Provenance]
    artifact_path: str | None = None
    rank: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class InteractionFingerprint:
    target_agi: str
    herbicide_selective_mutable: list[int]
    shared_protected: list[int]
    native_ligand_critical_protected: list[int]
    protected_by_context: list[int]
    second_shell_candidates: list[int]
    provenance: list[Provenance]


@dataclass(frozen=True)
class MutationCandidate:
    mutation: str
    agi: str
    structure_residue: int
    rationale: str
    classification: str
    provenance: list[Provenance]


@dataclass(frozen=True)
class ScorePacket:
    herbicide_escape_score: float
    native_ligand_retention_score: float
    functional_geometry_score: float
    fold_stability_score: float
    cofactor_or_complex_retention_score: float
    conservation_score: float
    pose_robustness_score: float
    method_agreement_score: float
    experimental_actionability_score: float
    provenance: list[Provenance]


@dataclass(frozen=True)
class EvaluationPacket:
    candidate: MutationCandidate
    scores: ScorePacket
    known_facts: list[str]
    predictions: list[str]
    assumptions: list[str]
    unresolved_uncertainty: list[str]
    recommendation: Literal["advance", "reject", "more_computation"]
    status: Literal[
        "DRAFT", "NEEDS_REVIEW", "APPROVED_FOR_ASSAY_PLANNING", "REJECTED", "NEEDS_MORE_COMPUTATION"
    ]
    provenance: list[Provenance]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class WorkflowRequest:
    target: TargetProtein
    herbicide: Ligand
    native_ligands: list[Ligand]
    functional_context: list[FunctionalPartner]
    protected_residues: set[int] = field(default_factory=set)


@dataclass(frozen=True)
class WorkflowResult:
    registry_entry: TargetRegistryEntry
    structures: list[StructureModel]
    poses: list[Pose]
    packets: list[EvaluationPacket]
    rejected_mutations: list[dict[str, str]]
    fingerprints: list[InteractionFingerprint]
