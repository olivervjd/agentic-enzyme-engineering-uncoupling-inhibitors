from __future__ import annotations

from ..agents.pipeline_agents import CandidateReviewAgent, ConstrainedMutationAgent, MultiOracleScoringAgent
from ..backends.interfaces import (
    AffinityPredictionBackend,
    ComplexModelingBackend,
    DockingBackend,
    RosalindReasoningBackend,
    StructurePredictionBackend,
)
from ..evals.deterministic_checks import validate_evaluation_packet
from ..registry.loader import TargetRegistry
from ..schemas.models import WorkflowRequest, WorkflowResult
from ..validators.input_validator import validate_request
from ..validators.provenance_validator import validate_provenance


class WorkflowOrchestrator:
    def __init__(
        self,
        registry: TargetRegistry,
        structure_backend: StructurePredictionBackend,
        docking_backend: DockingBackend,
        complex_backend: ComplexModelingBackend,
        affinity_backend: AffinityPredictionBackend,
        reasoning_backend: RosalindReasoningBackend,
    ) -> None:
        self.registry = registry
        self.structure_backend = structure_backend
        self.docking_backend = docking_backend
        self.complex_backend = complex_backend
        self.mutation_agent = ConstrainedMutationAgent()
        self.scoring_agent = MultiOracleScoringAgent(affinity_backend)
        self.review_agent = CandidateReviewAgent(reasoning_backend)

    def run(self, request: WorkflowRequest) -> WorkflowResult:
        entry = validate_request(request, self.registry)
        structures = self.structure_backend.predict_ensemble(request.target, entry)
        context_errors = [
            error for structure in structures for error in self.complex_backend.validate_context(structure, entry)
        ]
        if context_errors:
            raise ValueError(f"Structure context validation failed: {context_errors}")
        native_poses = [
            pose
            for ligand in request.native_ligands
            for pose in self.docking_backend.dock(structures, ligand)
        ]
        herbicide_poses = self.docking_backend.dock(structures, request.herbicide)
        candidates, rejected = self.mutation_agent.propose(request.target, request.protected_residues)
        facts = [
            f"Fixed registry pairing: {entry.agi} / {entry.herbicide}.",
            f"Mechanism class: {entry.mechanism_class}.",
            f"Required workflow: {entry.modeling_workflow}.",
        ]
        packets = [
            self.review_agent.review(
                candidate,
                self.scoring_agent.score(candidate, native_poses, herbicide_poses),
                facts,
            )
            for candidate in candidates
        ]
        validation_failures = [
            failure
            for packet in packets
            for failure in validate_evaluation_packet(packet, request.target.sequence, request.protected_residues)
        ]
        if validation_failures:
            raise ValueError(f"Evaluation packet validation failed: {validation_failures}")
        result = WorkflowResult(entry, structures, native_poses + herbicide_poses, packets, rejected)
        validate_provenance(result)
        return result
