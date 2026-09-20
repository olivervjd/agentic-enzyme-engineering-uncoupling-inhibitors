"""Mock subprocesses verify scheduling/provenance, without claiming Vina ran."""
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import threading
import time
import unittest
from herbicide_desensitization_agent.examples.epsps_vina import save
from herbicide_desensitization_agent.examples.epsps_vina_parallel import run_parallel
from herbicide_desensitization_agent.app.orchestrator.execution import ExecutionJournal,validate_journal

class ParallelVinaTests(unittest.TestCase):
    def test_bounded_parallel_merge_keeps_child_events_and_artifact_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp);campaign=base/'campaign';campaign.mkdir()
            save(campaign/'effective_model_inputs.json',{'records':[{'id':f'record{i}','mode':'structure','subject':'fixture','context':'herbicide','structure':'fixture.cif'} for i in range(4)]})
            active=0;peak=0;lock=threading.Lock();child_ids=[]
            def fake(command,**kwargs):
                nonlocal active,peak
                folder=Path(command[command.index('--output')+1]);identity=command[command.index('--record')+1]
                with lock:active+=1;peak=max(peak,active)
                folder.mkdir();job=folder/(identity+'--glyphosate');job.mkdir()
                artifact=job/'docking_results.json';save(artifact,{'fixture':'not molecular computation'})
                journal=ExecutionJournal(folder/'execution_journal.json')
                journal.execute('docking',lambda:{'fixture':True},output_artifacts=[artifact]);journal.finalize()
                child_ids.append(journal.data['stages'][0]['id'])
                save(folder/'campaign_results.json',{'records':[{'id':identity,'ligand':'glyphosate','execution_status':'COMPLETED','result':str(artifact)}]})
                time.sleep(.03)
                with lock:active-=1
                return SimpleNamespace(returncode=0)
            output=base/'output';result=run_parallel(campaign,output,workers=2,exclude_records=['record3'],runner=fake)
            self.assertEqual(peak,2);self.assertEqual(len(result['records']),3)
            merged=json.loads((output/'execution_journal.json').read_text());validate_journal(merged)
            self.assertEqual({r['source_event_id'] for r in merged['stages']},set(child_ids))
            self.assertEqual(len(merged['stages']),3)
            for row in result['records']:
                self.assertTrue(Path(row['result']).is_file())
                self.assertTrue((output/(row['id']+'--glyphosate')/'docking_results.json').is_file())
            for event in merged['stages']:
                self.assertIn('/_workers/',event['output_artifacts'][0]['path'])
                self.assertTrue(Path(event['output_artifacts'][0]['path']).is_file())

    def test_worker_failure_never_creates_completed_scientific_event(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp);campaign=base/'campaign';campaign.mkdir()
            save(campaign/'effective_model_inputs.json',{'records':[{'id':'fixture','mode':'structure','context':'native','structure':'fixture.cif'}]})
            result=run_parallel(campaign,base/'output',runner=lambda *a,**k:SimpleNamespace(returncode=9))
            self.assertEqual([r['execution_status'] for r in result['records']],['FAILED','FAILED'])
            journal=json.loads((base/'output'/'execution_journal.json').read_text())
            self.assertEqual(journal['stages'],[]);self.assertEqual(journal['status'],'FAILED')

if __name__=='__main__':unittest.main()
