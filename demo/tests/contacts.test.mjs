import test from 'node:test';
import assert from 'node:assert/strict';
import {contactScope,contactColumn} from '../src/contacts.js';
test('independent docking receptors and prediction votes are never pooled',()=>{
 const rows=[{method:'Boltz',ligand:'s3p',context:'native',contact_frequency:1},
 {method:'Vina',source_record:'receptor211',ligand:'s3p',contact_frequency:1},
 {method:'Vina',source_record:'receptor223',ligand:'s3p',contact_frequency:0}];
 const scope=contactScope(rows,'Vina','receptor223');
 assert.equal(scope.rows.length,1);
 assert.equal(scope.rows[0].contact_frequency,0);
 assert.deepEqual(scope.receptors,['receptor211','receptor223']);
});
test('same ligand in different chemical contexts has a separate matrix column',()=>{
 assert.notEqual(contactColumn({ligand:'s3p',context:'native'}),contactColumn({ligand:'s3p',context:'herbicide'}));
});
test('switching protein resets an unavailable receptor to an actual record',()=>{
 const scope=contactScope([{method:'Vina',source_record:'ecoli211'}],'Vina','arabidopsis211');
 assert.equal(scope.selectedReceptor,'ecoli211');assert.equal(scope.rows.length,1);
});
