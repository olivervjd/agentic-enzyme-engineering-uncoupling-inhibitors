"""Bounded independent Vina subprocesses; preserve actual child execution journals."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from .epsps_vina import jobs_from_manifest, digest, save
from ..app.orchestrator.execution import utc_now, validate_journal


def run_parallel(campaign_dir, output, *, workers=8, manifest_path=None, exclude_records=(), runner=subprocess.run):
    if isinstance(workers,bool) or not 1<=workers<=32: raise ValueError('Workers must be between 1 and 32')
    campaign_dir=Path(campaign_dir).resolve(); output=Path(output).resolve()
    manifest_path=Path(manifest_path).resolve() if manifest_path else campaign_dir/'effective_model_inputs.json'
    manifest=json.loads(manifest_path.read_text()); jobs=jobs_from_manifest(manifest)
    identities=list(dict.fromkeys(record['id'] for record,_ in jobs))
    unknown=set(exclude_records)-set(identities)
    if unknown: raise ValueError('Excluded record does not exist in eligible manifest: '+', '.join(sorted(unknown)))
    identities=[i for i in identities if i not in exclude_records]
    if not identities: raise ValueError('No eligible records remain')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',i) for i in identities): raise ValueError('Unsafe record identity')
    if output.exists(): raise FileExistsError('Fresh parallel output directory required')
    output.mkdir(parents=True); (output/'_workers').mkdir()
    started=utc_now(); monotonic=time.monotonic(); results={}; attempts=[]; journals=[]

    def worker(identity):
        folder=output/'_workers'/identity
        command=[sys.executable,'-m','herbicide_desensitization_agent.examples.epsps_vina',
                 '--campaign-dir',str(campaign_dir),'--manifest',str(manifest_path),
                 '--output',str(folder),'--record',identity]
        began=utc_now(); tick=time.monotonic()
        with (output/'_workers'/f'{identity}.log').open('w') as log:
            process=runner(command,stdout=log,stderr=subprocess.STDOUT,check=False)
        return identity,folder,{'record_id':identity,'command':command,'started_at':began,
            'finished_at':utc_now(),'wall_seconds':time.monotonic()-tick,'returncode':process.returncode,
            'scope':'Subprocess orchestration timing; excluded from scientific stage runtime totals'}

    def publish(finished=False):
        flat=[r for identity in identities for r in results.get(identity,[])]
        save(output/'campaign_results.json',{'records':flat,'manifest_sha256':digest(manifest_path),
            'prepare_only':False,'excluded_records':list(exclude_records),'pending_record_ids':[i for i in identities if i not in results],
            'orchestration':'Bounded independent subprocesses; one record per process, ligands serial within record'})
        events=[{**event,'id':journal['run_id']+':'+event['id'],'source_event_id':event['id'],'source_run_id':journal['run_id']} for _,journal in journals for event in journal['stages']]
        merged={'schema_version':'1.0','run_id':'parallel-'+output.name,'started_at':started,
            'finished_at':utc_now() if finished else None,
            'status':('FAILED' if any(r.get('execution_status')=='FAILED' for r in flat) else 'COMPLETED') if finished else 'RUNNING',
            'stages':events,'source_journals':[{'path':str(path),'sha256':digest(path)} for path,_ in journals],
            'merge_policy':'Original child events retained; merged ids qualified by source run id to avoid collisions, source_event_id preserved; timestamps, durations and artifact paths unchanged; nested seed journals remain artifacts to avoid double-counting ensemble duration'}
        validate_journal(merged); save(output/'execution_journal.json',merged)
        save(output/'parallel_orchestration.json',{'started_at':started,'finished_at':utc_now() if finished else None,
             'wall_seconds':time.monotonic()-monotonic,'workers':workers,'attempts':attempts,
             'scope':'Scheduling audit only; no fabricated scientific execution events'})

    publish()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,identity):identity for identity in identities}
        for future in as_completed(futures):
            identity,folder,attempt=future.result();attempts.append(attempt)
            result_path=folder/'campaign_results.json'
            records=json.loads(result_path.read_text())['records'] if result_path.is_file() else []
            for row in records:
                if row['id']!=identity: raise ValueError('Worker result identity mismatch')
                job_folder=folder/f"{identity}--{row['ligand']}"
                if job_folder.is_dir(): shutil.copytree(job_folder,output/job_folder.name)
            recorded={row['ligand'] for row in records}
            for record,target in jobs:
                if record['id']==identity and target not in recorded:
                    records.append({'id':identity,'ligand':target,'execution_status':'FAILED',
                        'reason':f"Worker exited with code {attempt['returncode']} before emitting this job result; no scientific completion is claimed"})
            if attempt['returncode'] and not any(r['execution_status']=='FAILED' for r in records):
                # Completed jobs stay completed; orchestration failure remains independently visible.
                attempt['warning']='Worker returned nonzero after writing its available completed jobs'
            results[identity]=records
            journal_path=folder/'execution_journal.json'
            if journal_path.is_file():
                journal=json.loads(journal_path.read_text());validate_journal(journal);journals.append((journal_path,journal))
            publish()
    publish(finished=True)
    return json.loads((output/'campaign_results.json').read_text())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-dir',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=8);parser.add_argument('--manifest',type=Path)
    parser.add_argument('--exclude-record',action='append',default=[])
    args=parser.parse_args(); result=run_parallel(args.campaign_dir,args.output,workers=args.workers,manifest_path=args.manifest,exclude_records=args.exclude_record)
    print(json.dumps({'jobs':len(result['records']),'failed':sum(r['execution_status']=='FAILED' for r in result['records'])}))
    if any(r['execution_status']=='FAILED' for r in result['records']):raise SystemExit(1)

if __name__=='__main__':main()
