from __future__ import annotations

from ..agents.pipeline_agents import (
    CandidateReviewAgent,
    ConstrainedMutationAgent,
    InteractionFingerprintAgent,
    MultiOracleScoringAgent,
)
from ..agents.mutation_scoring import pareto_rank
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
from ..storage.artifact_store import ArtifactStore
from ..validators.input_validator import validate_request
from ..validators.provenance_validator import validate_provenance
from ..validators.structure_context_validator import validate_structure_context
from ..visualization import MolstarArtifactRenderer


class WorkflowOrchestrator:
    def __init__(
        self,
        registry: TargetRegistry,
        structure_backend: StructurePredictionBackend,
        docking_backend: DockingBackend,
        complex_backend: ComplexModelingBackend,
        affinity_backend: AffinityPredictionBackend,
        reasoning_backend: RosalindReasoningBackend,
        artifact_store: ArtifactStore | None = None,
    ) -> None:
        self.registry = registry
        self.structure_backend = structure_backend
        self.docking_backend = docking_backend
        self.complex_backend = complex_backend
        self.mutation_agent = ConstrainedMutationAgent()
        self.fingerprint_agent = InteractionFingerprintAgent()
        self.scoring_agent = MultiOracleScoringAgent(affinity_backend)
        self.review_agent = CandidateReviewAgent(reasoning_backend)
        self.artifact_store = artifact_store

    def run(self, request: WorkflowRequest) -> WorkflowResult:
        entry = validate_request(request, self.registry)
        structures = self.structure_backend.predict_ensemble(request.target, entry)
        context_errors = [
            error
            for structure in structures
            for error in (
                self.complex_backend.validate_context(structure, entry)
                + validate_structure_context(structure, entry)
            )
        ]
        if context_errors:
            raise ValueError(f"Structure context validation failed: {context_errors}")
        native_poses = [
            pose
            for ligand in request.native_ligands
            for pose in self.docking_backend.dock(structures, ligand)
        ]
        herbicide_poses = self.docking_backend.dock(structures, request.herbicide)
        fingerprint = self.fingerprint_agent.compare(
            entry.agi, native_poses, herbicide_poses, request.protected_residues
        )
        candidates, rejected = self.mutation_agent.propose(request.target, request.protected_residues, fingerprint)
        facts = [
            f"Fixed registry pairing: {entry.agi} / {entry.herbicide}.",
            f"Mechanism class: {entry.mechanism_class}.",
            f"Required workflow: {entry.modeling_workflow}.",
        ]
        score_packets = [
            self.scoring_agent.score(candidate, native_poses, herbicide_poses, fingerprint)
            for candidate in candidates
        ]
        packets = [
            self.review_agent.review(
                candidate,
                scores,
                facts,
            )
            for candidate, scores in zip(candidates, score_packets)
        ]
        ranking = pareto_rank(candidates, score_packets)
        validation_failures = [
            failure
            for packet in packets
            for failure in validate_evaluation_packet(packet, request.target.sequence, request.protected_residues)
        ]
        if validation_failures:
            raise ValueError(f"Evaluation packet validation failed: {validation_failures}")
        result = WorkflowResult(
            entry, structures, native_poses + herbicide_poses, packets, rejected, [fingerprint], ranking
        )
        validate_provenance(result)
        if self.artifact_store:
            run_id = f"{entry.agi}-{entry.herbicide}".replace(",", "").replace(" ", "-")
            self.artifact_store.write_manifest(run_id, "structure_ensemble.json", structures)
            self.artifact_store.write_manifest(run_id, "pose_ensemble.json", native_poses + herbicide_poses)
            self.artifact_store.write_manifest(run_id, "interaction_fingerprint.json", fingerprint)
            self.artifact_store.write_manifest(run_id, "evaluation_packets.json", packets)
            self.artifact_store.write_manifest(run_id, "candidate_mutations.json", candidates)
            self.artifact_store.write_manifest(run_id, "rejected_mutations.json", rejected)
            self.artifact_store.write_manifest(run_id, "score_packets.json", score_packets)
            self.artifact_store.write_manifest(run_id, "pareto_ranking.json", ranking)
            renderer = MolstarArtifactRenderer(self.artifact_store)
            visualizations = [
                renderer.render(run_id, structure, index)
                for index, structure in enumerate(structures, start=1)
                if structure.artifact_path
            ]
            self.artifact_store.write_manifest(run_id, "visualizations.json", visualizations)
        return result
