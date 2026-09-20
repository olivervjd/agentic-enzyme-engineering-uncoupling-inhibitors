import unittest
from pathlib import Path
import tempfile
from herbicide_desensitization_agent.app.backends.calibration import sha
from herbicide_desensitization_agent.examples.export_fresh_epsps import (aggregate_contacts, binding_rows, validate_predictions,
    evaluate_readiness, nomination_readiness, function_diagnostic, baseline_availability, annotate_experimental_binding)


class FreshExportTests(unittest.TestCase):
    def test_crystal_binding_observation_preserves_separate_model_decision(self):
        bindings=[{'variant':subject,'ligand':ligand,'evidence_state':'Insufficient evidence','decision':'INSUFFICIENT_EVIDENCE'}
                  for subject in ('arabidopsis_wt','ecoli_wt','ecoli_G96A') for ligand in ('glyphosate','pep','s3p')]
        geometry=[{'subject':'ecoli_wt','reference':'1G6S','reference_sha256':'a'*64,'ligand_atom_mapping':{'glyphosate':{},'s3p':{}}},
                  {'subject':'ecoli_G96A','reference':'1MI4','reference_sha256':'b'*64,'ligand_atom_mapping':{'s3p':{}}}]
        result=annotate_experimental_binding(bindings,geometry)
        self.assertEqual(sum('experimental_evidence' in row for row in result),3)
        for row in result:
            self.assertEqual(row['decision'],'INSUFFICIENT_EVIDENCE')
            self.assertEqual(row['evidence_state'],'Insufficient evidence')
            if row['variant']=='arabidopsis_wt': self.assertNotIn('experimental_evidence',row)
            if 'experimental_evidence' in row:
                self.assertEqual(row['experimental_evidence']['evidence_ladder_state'],'Experimentally confirmed binding')
                self.assertIsNone(row['experimental_evidence']['Kd'])

    def test_available_baseline_is_not_reported_missing_but_acceptance_stays_uncalibrated(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); msa=root/'query.a3m'; msa.write_text('>query\nKDE\n>homolog\nKDA\n')
            manifest={'subjects':{'arabidopsis_wt':{'sequence':'KDE','offset':76,'msa':{'path':msa.name,'sha256':sha(msa)}}},
                      'seeds':[1,2,3],'registered_at':'2026-09-20T00:00:00+00:00',
                      'records':[{'mode':'structure','context':c,'structure':c+'.cif'} for c in ('apo','s3p_only','herbicide','native')]}
            baseline={'sources':{'at2022':{'url':'primary-at'},'ec2009':{'url':'primary-ec'},'catalytic':{'url':'primary-catalytic'}},
                      'mapping':{'residues':[{'source':1,'target':77}]},'protected_catalytic_positions':[{'canonical':77,'residue':'K','source':'catalytic'}],'controls':[]}
            comparisons=[{'subject':'arabidopsis_wt','context':c,'kind':'repeatability','query':str(i),'reference':str(j)} for c in ('apo','s3p_only','herbicide','native') for i,j in ((1,2),(1,3),(2,3))]
            comparisons.append({'subject':'arabidopsis_wt','context':'apo','kind':'experimental_accuracy_same_target_apo','query_tm_score':.8,'reference_sha256':'a'*64})
            docking=[{'source_record':'experimental','context':'herbicide','reference_scope':'experimental_redocking','independent_docking_seed_count':3,'provenance':{'sha256':'b'*64}}]
            readiness=evaluate_readiness(manifest,[],[],comparisons,baseline,root,docking)
            available=readiness['core_calibration']['specification']['baseline_evidence']
            for key in ('documented_msa','repeated_structure_prediction','experimental_structure_comparison','published_binding_site','catalytic_residues','biological_assembly','repeated_docking'):
                self.assertTrue(available[key]['validated'],key)
                self.assertFalse(available[key]['acceptance_calibrated'],key)
            self.assertFalse(any('WT baseline evidence incomplete' in reason for reason in readiness['errors']))
            self.assertFalse(readiness['eligible']); self.assertEqual(readiness['core_calibration']['thresholds'],{})
            msa.write_text('>query\nBAD\n')
            changed=baseline_availability(manifest,comparisons,baseline,root,docking)
            self.assertFalse(changed['documented_msa']['validated'])

    def test_real_core_gates_reject_unregistered_thresholds_without_fabricating_bounds(self):
        manifest={'seeds':[1,2,3],'registered_at':'2026-09-20T00:00:00+00:00','prospective_rules':{'contact_reproducibility_min':2/3}}
        readiness=evaluate_readiness(manifest,[],[],[],{'sources':{},'controls':[]},Path('/campaign'))
        self.assertFalse(readiness['eligible'])
        self.assertEqual(readiness['core_calibration']['thresholds'],{})
        self.assertIn('Preregistered metric specifications required',readiness['errors'])
        nomination=nomination_readiness([],readiness)
        self.assertFalse(nomination['candidate_generation_enabled'])
        self.assertIn('All preregistered differential-contact thresholds required',nomination['errors'])
        diagnostic=function_diagnostic([],readiness)
        self.assertEqual(diagnostic['decision'],'INSUFFICIENT_EVIDENCE')
        self.assertTrue(any(g['name']=='wt_calibration' and g['status']=='MISSING' for g in diagnostic['gates']))

    def test_native_s3p_uses_native_complex_and_excludes_affinity(self):
        manifest={'subjects':{'arabidopsis_wt':{'sequence':'A'*444,'offset':76}}}
        predictions=[{'id':f'{context}-{seed}','mode':'structure','subject':'arabidopsis_wt','context':context,'seed':seed}
                     for context in ('native','herbicide') for seed in (1,2,3)]
        predictions.append({'id':'affinity','mode':'affinity','subject':'arabidopsis_wt','context':'native','seed':1})
        def contact(ligand,context,seed):
            return {'variant':'arabidopsis_wt','residue':100,'ligand':ligand,'context':context,'pose_id':f'{context}-{seed}','min_distance':3.,'atom_pairs':[]}
        raw=[contact('s3p','native',seed) for seed in (1,2)] + [contact('pep','native',1),contact('glyphosate','herbicide',1)]
        # A separate affinity-mode record must not affect either native denominator.
        result=aggregate_contacts(raw,predictions,manifest,{})
        s3p=next(r for r in result if r['residue']==100 and r['ligand']=='s3p')
        self.assertEqual(s3p['replicates'],3); self.assertEqual(s3p['contact_frequency'],2/3)
        self.assertEqual(s3p['native_contact_frequencies'],{'pep':1/3,'s3p':2/3})
        catalytic=[r for r in result if r['residue']==407]
        self.assertEqual(len(catalytic),3)
        self.assertTrue(all(r['classification']=='CATALYTIC_OR_PROTECTED' for r in catalytic))

    def test_binding_uses_paired_affinity_seeds_and_no_confidence_interval(self):
        predictions=[]
        for subject,values in [('ecoli_wt',[1,3,6]),('ecoli_G96A',[2,2,5])]:
            for seed,value in enumerate(values):
                predictions.append({'mode':'affinity','id':subject+str(seed),'subject':subject,'ligand':'glyphosate','context':'herbicide','seed':seed,'predicted_pIC50':value,'affinity_sha256':'hash','ligands':[{'name':'glyphosate','model_smiles':'neutral'}]})
        predictions.append({'mode':'structure','predicted_pIC50':100})
        result=binding_rows(predictions)
        mutant=next(r for r in result if r['variant']=='ecoli_G96A')
        self.assertEqual([p['mutant_minus_wt'] for p in mutant['paired_seed_deltas']],[1,-1,-1])
        self.assertAlmostEqual(mutant['wt_relative_change'],-1/3)
        self.assertIsNone(mutant['interval']); self.assertEqual(mutant['replicates'],3)
        predictions[0]['ligands']=[{'model_smiles':'different'}]
        with self.assertRaisesRegex(ValueError,'chemistry'): binding_rows(predictions)

    def test_missing_and_duplicate_prediction_panel_fails(self):
        manifest={'subjects':dict.fromkeys(('arabidopsis_wt','ecoli_wt','ecoli_G96A')),'seeds':[1,2,3],'records':[]}
        for subject in manifest['subjects']:
            for seed in manifest['seeds']:
                for mode,contexts in [('structure',('apo','s3p_only','herbicide','native')),('affinity',('herbicide','native','native_s3p'))]:
                    for context in contexts:
                        manifest['records'].append({'id':f'{subject}-{seed}-{mode}-{context}','subject':subject,'seed':seed,'mode':mode,'context':context,'request':'request','request_sha256':'sha','ligands':[],'predicted_pIC50':1})
        validate_predictions(manifest,list(manifest['records']))
        with self.assertRaises(ValueError): validate_predictions(manifest,manifest['records'][:-1])
        with self.assertRaises(ValueError): validate_predictions(manifest,manifest['records']+[manifest['records'][0]])
        incomplete={**manifest,'records':manifest['records'][:-1]}
        with self.assertRaisesRegex(ValueError,'63-job'): validate_predictions(incomplete,incomplete['records'])
