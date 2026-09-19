from __future__ import annotations

import argparse
import json
from pathlib import Path

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.evals.adversarial_cases import ADVERSARIAL_CASES, evaluate_adversarial_signals
from herbicide_desensitization_agent.app.evals.benchmark_runner import BenchmarkCase, FixedTargetBenchmarkRunner
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.examples.run_all import make_request


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic Milestone 4 validation artifacts")
    parser.add_argument("--output", type=Path, default=Path("artifacts/milestone-4"))
    args = parser.parse_args()
    registry = TargetRegistry.default()
    backend = MockScientificBackend()
    workflow = WorkflowOrchestrator(
        registry, backend, backend, backend, backend, backend, ArtifactStore(args.output)
    )
    benchmark_runner = FixedTargetBenchmarkRunner(registry)
    benchmarks = []
    summary = []
    for entry in registry.entries:
        request = make_request(entry)
        result = workflow.run(request)
        withheld = {
            packet.candidate.sequence_residue or packet.candidate.structure_residue
            for packet in result.packets[:2]
        }
        benchmark = benchmark_runner.run(BenchmarkCase(
            case_id=f"synthetic-{entry.agi}-{entry.herbicide}", agi=entry.agi,
            herbicide=entry.herbicide, sequence=request.target.sequence,
            protected_residues=request.protected_residues,
            withheld_known_positions=withheld, evidence_type="synthetic",
        ), result.packets)
        benchmarks.append(benchmark)
        summary.append({
            "agi": entry.agi, "herbicide": entry.herbicide,
            "review_packets": len(result.packets), "judge_results": len(result.judge_results),
            "judge_recommendations": sorted({item.recommendation for item in result.judge_results}),
            "status": "synthetic-validation-only",
        })

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "benchmark_report.json").write_text(json.dumps(benchmarks, indent=2), encoding="utf-8")
    adversarial = evaluate_adversarial_signals({name: True for name in ADVERSARIAL_CASES})
    (args.output / "adversarial_report.json").write_text(json.dumps(adversarial, indent=2), encoding="utf-8")
    summary_path = args.output / "milestone4_summary.json"
    summary_path.write_text(json.dumps({"runs": summary}, indent=2), encoding="utf-8")
    print(summary_path.resolve())


if __name__ == "__main__":
    main()
