"""Execute evidence agents and assemble current-campaign journals without invented work."""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

from ..app.backends.calibration import sha, write
from ..app.backends.credentials import resolve_api_key
from ..app.backends.openai_json import DEFAULT_AGENT_MODEL, OpenAIJSONTransport
from ..app.agents.evidence import CitedEvidenceAgent
from ..app.registry.loader import TargetRegistry
from ..app.orchestrator.execution import ExecutionJournal, validate_journal, utc_now


def literature(root, registry_path, source_dir):
    root, registry_path, source_dir = map(Path, (root, registry_path, source_dir))
    output = root / 'evidence'
    output.mkdir(parents=True, exist_ok=True)
    journal = ExecutionJournal(output / 'execution_journal.json')
    baseline = json.loads(registry_path.read_text())
    entry = TargetRegistry.default().get('AT2G45300', 'glyphosate')
    sources = []
    for key, source in baseline['sources'].items():
        # The artifact is explicitly a curation snapshot, not a downloaded paper.
        sources.append({'id': key, 'url': source['url'], 'scope': 'Source-linked curated baseline; applicability and missing endpoints retained',
                        'artifact': 'curated_baseline.json', 'sha256': sha(registry_path),
                        'title': source.get('title', source.get('finding', key)),
                        'artifact_type': 'curated_source_record', 'curated_record': source})
    (output / 'curated_baseline.json').write_bytes(registry_path.read_bytes())
    raw_sources = []
    for filename in ('PMC8983318.xml', 'Funke09.pdf', '7PXY.cif', 'chemical_states_ph74.json'):
        path = source_dir / filename
        # Missing primary files fail the journaled literature inputs, while an
        # available chemical-state input can still be validated independently.
        raw_sources.append(path)
    claims = [{'text': 'Curated EPSPS baseline (all numerical quantities originate from cited primary records, not this synthesis): ' + json.dumps(baseline, allow_nan=False),
               'source_ids': list(baseline['sources'])}]
    evidence = {'target_agi': entry.agi, 'sources': sources, 'claims': claims,
                'limitations': ['Homolog controls do not establish the same phenotype in Arabidopsis',
                    'Published kinetics remain Ki, Km, kcat or condition-dependent IC50; none are Kd',
                    'Ligand protonation hypotheses are not measured bound-state populations']}
    write(output / 'evidence_inputs.json', evidence)
    def synthesize():
        api_key, credential_source = resolve_api_key(use_codex=True)
        transport = OpenAIJSONTransport(DEFAULT_AGENT_MODEL, api_key=api_key)
        agent = CitedEvidenceAgent(output / 'evidence_inputs.json', transport)
        result = agent.synthesize(entry)
        result.update(model=DEFAULT_AGENT_MODEL, credential_source=credential_source,
                      scientific_values_source='Immutable curated baseline and primary artifacts',
                      primary_artifacts=[{'filename': p.name, 'sha256': sha(p)} for p in raw_sources])
        write(output / 'literature_synthesis.json', result)
        return result
    def chemical_validation():
        from rdkit import Chem
        document = json.loads((source_dir / 'chemical_states_ph74.json').read_text())
        rows = document['ligands']
        states = {}
        for name, row in rows.items():
            smiles = row['selected_campaign_state']['smiles']
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or any(flag == '?' for _, flag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)):
                raise ValueError('Chemical state failed identity/stereo check')
            states[name] = {'smiles': Chem.MolToSmiles(mol), 'formal_charge': Chem.GetFormalCharge(mol)}
        result = {'states': states, 'status': 'IDENTITY_VALIDATED_MICROSTATE_ASSUMED',
                  'ph': 7.4, 'microstate_experimentally_validated': False,
                  'uncertainty': 'Enumerated alternatives retained; populations and bound-state pKa unmeasured'}
        write(output / 'chemical_validation.json', result)
        return result
    failures = []
    try:
        try:
            journal.execute('literature', synthesize, inputs={'target': entry.agi, 'model': DEFAULT_AGENT_MODEL},
                input_artifacts=[registry_path, output / 'evidence_inputs.json', *raw_sources],
                output_artifacts=[output / 'literature_synthesis.json'],
                tools=[{'name': 'OpenAI Responses API', 'requested_model': DEFAULT_AGENT_MODEL, 'provider_revision': 'not returned by transport'}])
        except Exception as exc:
            failures.append({'stage': 'literature', 'error_type': type(exc).__name__})
        try:
            journal.execute('chemical_state', chemical_validation, inputs={'ph_assumption': 7.4},
                input_artifacts=[source_dir / 'chemical_states_ph74.json'], output_artifacts=[output / 'chemical_validation.json'])
        except Exception as exc:
            failures.append({'stage': 'chemical_state', 'error_type': type(exc).__name__})
        if failures:
            write(output / 'execution_failures.json', failures)
            raise RuntimeError('Evidence-stage failure; independent literature and chemistry stages were both attempted')
    finally:
        journal.finalize()


def merge_journals(paths, destination):
    """Copy original executed events; do not turn imports/skips into execution."""
    sources, events, started, finished = [], [], [], []
    for number, path in enumerate(map(Path, paths), 1):
        data = json.loads(path.read_text())
        validate_journal(data)
        started.append(data['started_at'])
        finished.append(data.get('finished_at'))
        sources.append({'path': str(path), 'sha256': sha(path), 'run_id': data['run_id']})
        event_ids = {row['id'] for row in data['stages']}
        by_stage = {}
        for row in data['stages']:
            by_stage.setdefault(row['stage'], []).append(row['id'])
        for event in data['stages']:
            copy = deepcopy(event)
            copy['source_event_id'] = event['id']
            copy['source_run_id'] = data['run_id']
            copy['source_stage'] = event['stage']
            if copy['stage'] == 'stability_scoring':
                copy['stage'] = 'stability'
            copy['id'] = f'{number}:{event["id"]}'
            copy['source_dependencies'] = list(event.get('dependencies', []))
            resolved, unresolved = [], []
            for dependency in copy['source_dependencies']:
                if dependency in event_ids:
                    resolved.append(f'{number}:{dependency}')
                elif dependency in by_stage:
                    resolved.extend(f'{number}:{identity}' for identity in by_stage[dependency])
                else:
                    unresolved.append({'source_run_id': data['run_id'], 'reference': dependency})
            copy['dependencies'] = list(dict.fromkeys(resolved))
            if unresolved:
                copy['unresolved_source_dependencies'] = unresolved
            events.append(copy)
    if not sources:
        raise ValueError('At least one actual source journal is required')
    status = ('RUNNING' if any(x is None for x in finished) else
              'COMPLETED_WITH_FAILURES' if any(x['execution_status'] == 'FAILED' for x in events) else
              'COMPLETED_WITH_SKIPS' if any(x['execution_status'].startswith('SKIPPED') for x in events) else 'COMPLETED')
    data = {'schema_version': '1.0', 'run_id': 'epsps-fresh-20260920', 'started_at': min(started),
            'finished_at': max(finished) if all(finished) else None, 'status': status, 'stages': events,
            'sources': sources, 'assembled_at': utc_now(), 'assembly_is_not_new_scientific_execution': True}
    validate_journal(data)
    write(destination, data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['literature', 'merge'])
    parser.add_argument('--root', type=Path)
    parser.add_argument('--registry', type=Path)
    parser.add_argument('--source-dir', type=Path)
    parser.add_argument('--journals', nargs='+', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.action == 'literature':
        literature(args.root, args.registry, args.source_dir)
    else:
        merge_journals(args.journals, args.output)


if __name__ == '__main__':
    main()
