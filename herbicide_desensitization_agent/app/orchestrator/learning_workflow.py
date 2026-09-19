from __future__ import annotations

from ..agents.learning_agent import LearningAgent
from ..schemas.models import AssayResult, WorkflowRequest, WorkflowResult
from ..storage.artifact_store import ArtifactStore


class LearningWorkflow:
    def __init__(self, agent: LearningAgent, artifact_store: ArtifactStore) -> None:
        self.agent = agent
        self.artifact_store = artifact_store

    def run(
        self,
        request: WorkflowRequest,
        result: WorkflowResult,
        assays: list[AssayResult],
        limit: int = 3,
    ):
        entry = result.registry_entry
        if request.target.agi != entry.agi or request.herbicide.name.casefold() != entry.herbicide.casefold():
            raise ValueError("Learning request does not match workflow result")
        report = self.agent.recalibrate(
            entry.agi, entry.herbicide, result.packets, assays, request.target.sequence
        )
        next_round = self.agent.select_next_round(result.packets, assays, report, limit)
        run_id = f"{entry.agi}-{entry.herbicide}".replace(",", "").replace(" ", "-")
        self.artifact_store.write_manifest(run_id, "assay_results.json", assays)
        self.artifact_store.write_manifest(run_id, "model_recalibration_report.json", report)
        self.artifact_store.write_manifest(run_id, "next_round_candidates.json", next_round)
        return report, next_round
