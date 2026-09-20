"""Exporter statuses must prove execution, never infer it from existing results."""
import unittest
import tempfile
from pathlib import Path

from herbicide_desensitization_agent.app.orchestrator.execution import ExecutionJournal
from herbicide_desensitization_agent.examples.execution_table import execution_table


def row_for(rows, stage):
    return next(row for row in rows if row['agent_id'] == stage)


class ExecutionTableTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.tmp_path = Path(self.tmp.name)

    def test_no_journal_does_not_invent_completed_agents(self):
        rows, summary = execution_table(imported_sources=[{'path':'recorded-predictions.json','origin':'recorded_run'}])
        assert len(rows) == 14
        assert {r['execution_status'] for r in rows} == {'not_started'}
        assert {r['evidence_status'] for r in rows} == {'NOT_ASSESSED'}
        assert all(r['runtime'] is None and r['outputs'] == [] and r['tools_models'] == [] for r in rows)
        assert summary['journal'] is None
        assert summary['counts'] == {'not_started':14}


    def test_real_job_success_does_not_mean_scientific_validation(self):
        tmp_path = self.tmp_path
        journal = ExecutionJournal(tmp_path/'journal.json', run_id='test-run')
        result = tmp_path/'result.json'
        def work():
            result.write_text('{"computed":true}')
            return {'result':'result.json'}
        journal.execute('calibration', work, output_artifacts=[result], evidence_status='INSUFFICIENT_EVIDENCE')
        journal.finalize()
        rows, summary = execution_table(journal.data, journal_path=journal.path)
        row = row_for(rows, 'calibration')
        assert row['execution_status'] == 'succeeded'
        assert row['evidence_status'] == 'INSUFFICIENT_EVIDENCE'
        assert row['runtime'] == journal.stages[0]['runtime_seconds']
        assert row['outputs'][0]['sha256'] == journal.stages[0]['output_artifacts'][0]['sha256']
        assert row['provenance']['event_ids'] == [journal.stages[0]['id']]
        assert summary['run_id'] == 'test-run'
        assert len(summary['journal']['sha256']) == 64


    def test_verified_import_is_not_new_structure_execution(self):
        tmp_path = self.tmp_path
        prior = tmp_path/'prior.cif'
        prior.write_text('Recorded structure fixture')
        journal = ExecutionJournal()
        journal.import_artifacts('structure',[prior], evidence_status='AVAILABLE')
        journal.finalize()
        rows, _ = execution_table(journal.data)
        row = row_for(rows,'structure')
        assert row['execution_status'] == 'not_started'
        assert row['execution_mode'] == 'imported'
        assert row['evidence_status'] == 'AVAILABLE'
        assert row['runtime'] is None
        assert row['journal_execution_status'] == ['IMPORTED']


    def test_skipped_stages_and_failed_jobs_are_not_completed(self):
        journal = ExecutionJournal()
        journal.skip('mutation_design','SKIPPED_DEPENDENCY','Calibration was insufficient')
        journal.skip('stability','SKIPPED_NO_CANDIDATES','No eligible candidate')
        with self.assertRaises(ValueError):
            journal.execute('docking', lambda: (_ for _ in ()).throw(ValueError('unit test failure')))
        journal.finalize()
        rows, summary = execution_table(journal.data)
        assert row_for(rows,'mutation_design')['execution_status'] == 'skipped'
        assert row_for(rows,'stability')['runtime'] is None
        assert row_for(rows,'docking')['execution_status'] == 'failed'
        assert summary['counts'].get('succeeded', 0) == 0


    def test_partial_stage_does_not_claim_complete_execution(self):
        journal = ExecutionJournal()
        journal.execute('affinity', lambda: {'result':'first replicate'})
        journal.skip('affinity','SKIPPED_DEPENDENCY','Second replicate input missing')
        journal.finalize()
        rows, _ = execution_table(journal.data)
        assert row_for(rows,'affinity')['execution_status'] == 'skipped'
        assert len(row_for(rows,'affinity')['execution_records']) == 2


    def test_running_job_has_no_finished_runtime(self):
        journal = ExecutionJournal()
        journal.start('structure')
        rows, _ = execution_table(journal.data)
        assert row_for(rows,'structure')['execution_status'] == 'running'
        assert row_for(rows,'structure')['runtime'] is None


    def test_forged_completed_status_without_timing_is_rejected(self):
        journal = ExecutionJournal()
        event = journal.start('structure')
        event.update(execution_status='COMPLETED',outputs_sha256='0'*64)
        with self.assertRaises(ValueError):
            execution_table(journal.data)


    def test_unknown_stage_not_silently_assigned_to_agent(self):
        journal = ExecutionJournal()
        journal.execute('unmapped_preflight', lambda: {'available':False})
        journal.finalize()
        rows, summary = execution_table(journal.data)
        assert all(r['execution_status'] == 'not_started' for r in rows)
        assert summary['unmapped_events'][0]['stage'] == 'unmapped_preflight'

if __name__ == "__main__":
    unittest.main()
