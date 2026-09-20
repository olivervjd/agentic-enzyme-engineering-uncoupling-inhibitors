// Execution success is never inferred from a scientific evidence verdict.
const states = new Set(['not_started','running','succeeded','failed','skipped']);
export function executionState(record) {
  if (!record || !states.has(record.execution_status)) return 'not_started';
  if (record.execution_mode === 'imported') return 'not_started';
  return record.execution_status;
}
export function evidenceState(record) {
  return record?.evidence_status || 'NOT_ASSESSED';
}
export function executionCounts(records) {
  const counts={not_started:0,running:0,succeeded:0,failed:0,skipped:0};
  for(const record of records)counts[executionState(record)]++;
  return counts;
}
export function executionLabel(record) {
  return record?.execution_mode==='imported'?'prior evidence imported':executionState(record).replaceAll('_',' ');
}
