from __future__ import annotations

import argparse
import json
from pathlib import Path

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic Milestone 3 validation artifacts")
    parser.add_argument("--output", type=Path, default=Path("artifacts/milestone-3"))
    args = parser.parse_args()
    registry = TargetRegistry.default()
    backend = MockScientificBackend()
    workflow = WorkflowOrchestrator(
        registry, backend, backend, backend, backend, backend, ArtifactStore(args.output)
    )
    summary = []
    for entry in registry.entries:
        result = workflow.run(make_request(entry))
        summary.append({
            "agi": entry.agi,
            "herbicide": entry.herbicide,
            "candidate_count": len(result.packets),
            "pareto_front_1": [item.mutation for item in result.pareto_ranking if item.front == 1],
            "status": "synthetic-validation-only",
        })
    path = args.output / "milestone3_summary.json"
    path.write_text(json.dumps({"runs": summary}, indent=2), encoding="utf-8")
    print(path.resolve())


if __name__ == "__main__":
    main()
