from __future__ import annotations

import argparse
import json
from pathlib import Path

from herbicide_desensitization_agent.app.agents.learning_agent import LearningAgent
from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.orchestrator.learning_workflow import LearningWorkflow
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import AssayResult, Provenance
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


SYNTHETIC = [Provenance(
    "example://milestone-5-synthetic-assay", "deterministic-fixture", "synthetic",
    "Not experimental evidence",
)]


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic Milestone 5 learning-loop artifacts")
    parser.add_argument("--output", type=Path, default=Path("artifacts/milestone-5"))
    args = parser.parse_args()
    registry = TargetRegistry.default()
    backend = MockScientificBackend()
    store = ArtifactStore(args.output)
    workflow = WorkflowOrchestrator(registry, backend, backend, backend, backend, backend)
    learning = LearningWorkflow(LearningAgent(registry), store)
    summary = []
    for entry in registry.entries:
        request = make_request(entry)
        result = workflow.run(request)
        mutation = result.packets[0].candidate.mutation
        measurements = {
            "herbicide_response_fold_change": 2.0,
            "protein_expression_fraction": 0.85,
            "thermal_stability_delta_c": -0.5,
        }
        if entry.agi == "ATCG00020":
            measurements["photosynthetic_performance_fraction"] = 0.8
        elif entry.agi == "AT3G62980":
            measurements["native_signaling_fraction"] = 0.8
        else:
            measurements["native_activity_fraction"] = 0.8
        units = {name: ("degC" if name == "thermal_stability_delta_c" else "fold" if name.endswith("fold_change") else "fraction") for name in measurements}
        assay = AssayResult(
            f"synthetic-{entry.agi}-1", entry.agi, entry.herbicide, mutation,
            measurements, units, SYNTHETIC, "Synthetic learning-loop validation record",
        )
        report, next_round = learning.run(request, result, [assay])
        summary.append({
            "agi": entry.agi, "herbicide": entry.herbicide, "assayed_mutation": mutation,
            "next_round": [candidate.mutation for candidate in next_round],
            "assay_count": report.assay_count, "status": "synthetic-validation-only",
        })
    path = args.output / "milestone5_summary.json"
    path.write_text(json.dumps({"runs": summary}, indent=2), encoding="utf-8")
    print(path.resolve())


if __name__ == "__main__":
    main()
