"""Combine completed docking artifacts while preserving each originating campaign."""
import argparse
import json
from pathlib import Path
import shutil

from ..app.backends.calibration import sha, write
from .fresh_external_evidence import load_docking


def combine(sources, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError('Use a new combined-artifact directory')
    output.mkdir(parents=True)
    records, provenance, identities = [], [], set()
    for source in map(Path, sources):
        original = json.loads(source.read_text())
        if original.get('pending_record_ids'):
            raise ValueError('A source docking campaign is still running')
        jobs, failures, artifacts = load_docking(source)
        if failures:
            raise ValueError('Cannot describe failed source jobs as completed docking')
        provenance.append({'path': str(source), 'sha256': sha(source), 'manifest_sha256': original.get('manifest_sha256')})
        for path, job in jobs:
            record = job['preparation']['record_id']
            ligand = job['preparation']['docked_ligand']
            identity = (record, ligand)
            if identity in identities:
                raise ValueError('Duplicate docking identity across source campaigns')
            identities.add(identity)
            destination = output / f'{record}--{ligand}'
            destination.mkdir()
            for artifact in artifacts:
                if artifact.parent == path.parent:
                    shutil.copy2(artifact, destination / artifact.name)
            records.append({'id': record, 'ligand': ligand, 'execution_status': 'COMPLETED',
                            'result': str(destination / 'docking_results.json'), 'poses': len(job['rows']),
                            'originating_result_sha256': sha(path), 'originating_campaign_sha256': sha(source)})
    result = {'records': records, 'source_campaigns': provenance, 'pending_record_ids': [],
              'prepare_only': False, 'assembly_scope': 'Copies existing current-campaign artifacts; no new docking invocation claimed'}
    write(output / 'campaign_results.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = combine(args.sources, args.output)
    print(json.dumps({'combined_jobs': len(result['records'])}))
