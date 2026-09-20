"""Review, package and journal the executed EPSPS campaign; never invent assays."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from ..app.agents.learning_agent import LearningAgent
from ..app.backends.calibration import sha, write
from ..app.backends.credentials import resolve_api_key
from ..app.orchestrator.execution import ExecutionJournal
from ..app.registry.loader import TargetRegistry
from .epsps_evidence_stages import merge_journals
from .execution_table import execution_table
from .export_evidence_dashboard import table_export
from .package_evidence_dashboard import package
from .review_evidence_bundle import run_review


def finalize(root, bundle, dist, journals):
    root, bundle, dist = map(Path, (root, bundle, dist))
    journal = ExecutionJournal(root / 'publication_execution_journal.json')
    data_dir = bundle / 'data'
    manifest = json.loads((root / 'effective_model_inputs.json').read_text())
    evidence_path = data_dir / 'evidence_system.json'

    def intake():
        # Execute input validation. Empty input cannot trigger recalibration.
        LearningAgent(TargetRegistry.default()).validate_assays([], [], manifest['sequence'])
        result = {'status': 'AWAITING_ASSAYS', 'assay_count': 0, 'model_updated': False,
            'intake_validation_executed': True, 'reason': 'No new measured candidate assay records were supplied',
            'published_control_kinetics_are_not_new_candidate_assays': True}
        write(data_dir / 'learning_readiness.json', result)
        return result
    journal.execute('learning_intake', intake, inputs={'assays': [], 'reviewed_candidates': []},
        output_artifacts=[data_dir / 'learning_readiness.json'], evidence_status='INSUFFICIENT_EVIDENCE')
    journal.skip('learning_loop', 'SKIPPED_NO_INPUT', 'No new measured candidate assays; model fitting cannot execute',
                 inputs={'assay_count': 0}, dependencies=['learning_intake'])
    key, credential_source = resolve_api_key(use_codex=True)
    snapshot = bundle / 'raw/prepublication_evidence.json'
    shutil.copy2(evidence_path, snapshot)

    def review():
        source = json.loads(snapshot.read_text())
        result = run_review(source, api_key=key, credential_source=credential_source,
            source_provenance={'path': str(snapshot), 'sha256': sha(snapshot)})
        write(data_dir / 'prepublication_model_review.json', result)
        if result['status'] != 'COMPLETED':
            raise RuntimeError('Model review or judge did not execute successfully')
        return result
    journal.execute('evidence_review', review, inputs={'model': 'gpt-5.6-luna', 'source': str(snapshot)},
        input_artifacts=[snapshot], output_artifacts=[data_dir / 'prepublication_model_review.json'])

    def build():
        package(bundle, dist)
        return {'offline_report_created': True, 'figures': [p.name for p in (bundle / 'figures').glob('*')]}
    # Figure byte hashes remain stable when the final execution table is assembled.
    journal.execute('dashboard', build, inputs={'bundle': str(bundle), 'dist': str(dist)},
        input_artifacts=[snapshot, data_dir / 'structures.json', dist / 'index.html'],
        output_artifacts=lambda _: list((bundle / 'figures').glob('*')))
    journal.finalize()
    master_path = bundle / 'raw/execution_journal.json'
    combined = merge_journals([*journals, bundle / 'analysis_journal.json', journal.path], master_path)
    agents, execution = execution_table(combined, journal_path=master_path)
    if any(row['execution_status'] in {'failed', 'running', 'not_started'} for row in agents):
        raise RuntimeError('Required stages failed or have no recorded completed execution/explicit skip; do not publish completion')
    evidence = json.loads(evidence_path.read_text())
    evidence['tables']['agent_execution'] = agents
    evidence['execution'] = execution
    evidence['learning_readiness'] = json.loads((data_dir / 'learning_readiness.json').read_text())
    write(evidence_path, evidence)
    table_export(data_dir, 'agent_execution', agents, evidence['table_legends']['agent_execution'])
    workflow_manifest = json.loads((bundle / 'workflow-manifest.json').read_text())
    workflow_manifest.update(execution=execution, nodes=agents, current_campaign_complete=True,
        completion_scope='All available computational stages executed; gated mutation generation and assay-dependent learning retain explicit skips')
    write(bundle / 'workflow-manifest.json', workflow_manifest)
    # Execution metadata changed, so review the final immutable snapshot anew.
    # The first review remains journaled with its exact preserved source.
    result = run_review(evidence, api_key=key, credential_source=credential_source,
        source_provenance={'path': str(evidence_path), 'sha256': sha(evidence_path)})
    write(data_dir / 'evidence_model_review.json', result)
    if result['status'] != 'COMPLETED':
        raise RuntimeError('Final model commentary did not complete')
    package(bundle, dist, render_figures=False)
    return {'execution': execution['counts'], 'scientific_status': evidence['status'],
            'report': str(bundle / 'report.html')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--dist', type=Path, required=True)
    parser.add_argument('--journals', nargs='+', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(finalize(args.root, args.bundle, args.dist, args.journals)))


if __name__ == '__main__':
    main()
