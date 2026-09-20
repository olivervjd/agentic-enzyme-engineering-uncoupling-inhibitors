import test from 'node:test';
import assert from 'node:assert/strict';
import {executionState,evidenceState,executionCounts,executionLabel} from '../src/execution.js';
test('legacy completed and validated evidence cannot imply execution',()=>{
  assert.equal(executionState({status:'COMPLETED'}),'not_started');
  assert.equal(executionState({evidence_status:'VALIDATED'}),'not_started');
  assert.equal(executionState({execution_status:'succeeded',execution_mode:'imported'}),'not_started');
});
test('actual succeeded job retains insufficient scientific evidence',()=>{
  const row={execution_status:'succeeded',execution_mode:'executed',evidence_status:'INSUFFICIENT_EVIDENCE'};
  assert.equal(executionState(row),'succeeded');
  assert.equal(evidenceState(row),'INSUFFICIENT_EVIDENCE');
});
test('skips, failures, running jobs and imports are counted independently',()=>{
  assert.deepEqual(executionCounts([{execution_status:'skipped'},{execution_status:'failed'},{execution_status:'running'},{execution_status:'succeeded',execution_mode:'executed'},{execution_status:'succeeded',execution_mode:'imported'},{status:'COMPLETED'}]),{not_started:2,running:1,succeeded:1,failed:1,skipped:1});
  assert.equal(executionLabel({execution_mode:'imported'}),'prior evidence imported');
});
