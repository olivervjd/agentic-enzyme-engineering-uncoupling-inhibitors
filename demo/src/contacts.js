export const contactMethod = row => row.method || (row.model_mode === 'structure' ? 'Structure-only complex predictions' : '') || (row.methods_supporting || []).join(' + ') || (row.model_mode === 'docking' ? 'Docking (method unrecorded)' : 'Recorded complex predictions');
export const contactColumn = row => JSON.stringify([String(row.ligand || row.ligand_id || '').toLowerCase(), row.context || 'Context unrecorded', row.source_record || '']);
export function contactScope(rows, method, receptor) {
  const methods = [...new Set(rows.map(contactMethod))];
  const selectedMethod = methods.includes(method) ? method : methods[0];
  const matching = rows.filter(row => contactMethod(row) === selectedMethod);
  const receptors = [...new Set(matching.map(row => row.source_record).filter(Boolean))];
  const selectedReceptor = receptors.includes(receptor) ? receptor : receptors[0];
  return {methods, selectedMethod, receptors, selectedReceptor,
    rows: matching.filter(row => !selectedReceptor || row.source_record === selectedReceptor)};
}
