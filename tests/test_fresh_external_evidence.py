import json
from pathlib import Path
import tempfile
import unittest

from herbicide_desensitization_agent.examples.fresh_external_evidence import docking_tables, load_docking, stability_tables
from herbicide_desensitization_agent.examples.export_fresh_epsps import evaluate_readiness


class ExternalEvidenceTests(unittest.TestCase):
    def docking_fixture(self, root):
        jobdir=root/'wt-herbicide--glyphosate'; jobdir.mkdir()
        poses=[{'record_id':'wt-herbicide','subject':'arabidopsis_wt','context':'herbicide','ligand':'glyphosate',
                'seed':seed,'pose_rank':rank,'pose_id':f'{seed}:{rank}','vina_score_kcal_mol':-4.+rank/10,
                'reference_pose_rmsd_angstrom':float(rank),'reference_scope':'structure'} for seed in (211,223,227) for rank in (1,2,3)]
        contacts=[{'variant':'arabidopsis_wt','context':'herbicide','ligand':'glyphosate','chain':'A','structure_residue':24,
                   'residue':100,'replicate_count':3,'pose_count':9,'seeds_supporting_top_pose':[211,223],
                   'poses_supporting':['211:1','223:1'],'contact_frequency':2/3,'uncertainty':'Observed sampling fractions'}]
        job={'status':'COMPLETED','rows':poses,'pose_clusters':[{'cluster':0,'independent_seeds':[211,223,227],'fraction_seeds':1}],
             'clustering_threshold_angstrom':2.,'limitations':['Not thermodynamic affinity'],
             'preparation':{'record_id':'wt-herbicide','docked_ligand':'glyphosate'},'residue_contact_frequencies':contacts,'pose_contacts':[]}
        (jobdir/'docking_results.json').write_text(json.dumps(job))
        campaign=root/'campaign_results.json'
        campaign.write_text(json.dumps({'records':[{'id':'wt-herbicide','ligand':'glyphosate','execution_status':'COMPLETED',
                                                    'result':'/remote/run/wt-herbicide--glyphosate/docking_results.json'}]}))
        return campaign

    def test_pose_counts_do_not_become_independent_replicates_or_affinity_agreement(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.docking_fixture(Path(folder)); bindings=[{'variant':'arabidopsis_wt','ligand':'glyphosate','estimate':2.}]
            result,artifacts=docking_tables(path,bindings)
            summary=bindings[0]['docking_evidence'][0]
            self.assertEqual(summary['returned_pose_count'],9)
            self.assertEqual(summary['independent_docking_seed_count'],3)
            self.assertFalse(summary['chemical_context_equivalent_to_boltz_affinity'])
            self.assertEqual(bindings[0]['estimate'],2.)
            self.assertEqual(result['contacts'][0]['replicates'],3)
            self.assertEqual(result['contacts'][0]['contact_frequency'],2/3)
            self.assertEqual(result['controls'][0]['pass_fail'],'INSUFFICIENT_EVIDENCE')
            self.assertIn('not a redocking',summary['cluster_cutoff_scope'])

    def test_campaign_result_cannot_escape_job_allowlist(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.docking_fixture(Path(folder)); data=json.loads(path.read_text())
            data['records'][0]['result']='/outside/credentials.json'; path.write_text(json.dumps(data))
            with self.assertRaises(ValueError): load_docking(path)

    def test_experimental_redocking_scope_and_fallback_evidence_ladder(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.docking_fixture(Path(folder))
            source=path.parent/'wt-herbicide--glyphosate/docking_results.json'
            data=json.loads(source.read_text())
            for row in data['rows']: row['reference_scope']='experimental_redocking'
            source.write_text(json.dumps(data)); bindings=[]
            result,_=docking_tables(path,bindings)
            self.assertEqual(result['controls'][0]['experimental_agreement'],'Experimental pose comparison')
            self.assertEqual(bindings[0]['evidence_state'],'Insufficient evidence')
            self.assertIn('Docked poses',bindings[0]['pose_hypothesis'])

    def test_stability_is_diagnostic_and_actual_core_gates_block_nomination(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'stability_results.json'
            path.write_text(json.dumps({'method':'DDGun3D','unit':'kcal/mol','protocol':'custom_ColabFoldMSA_DDGun3D',
                'mutations':[{'mutation':'G177A','canonical_position':177,'scope':'known_control','folding_ddg_kcal_mol':1.2,
                    'ddgun3d_unfolding_ddg_kcal_mol':-1.2,'uncertainty_interval':None,'replicates':1,'calibration_status':'UNCALIBRATED_CUSTOM_MSA_PROTOCOL'}],
                'limitations':['Custom profile not calibrated'],'output_artifacts':[]}))
            readiness=evaluate_readiness({'seeds':[1,2,3],'registered_at':'2026-09-20T00:00:00+00:00'},[],[],[],{'sources':{},'controls':[]},Path(folder))
            result,_=stability_tables(path,readiness,[])
            row=result['mutations'][0]
            self.assertEqual(row['ddg'],1.2); self.assertFalse(row['new_nomination'])
            self.assertIsNone(row['uncertainty_interval']); self.assertIsNone(row['structure_retention'])
            self.assertEqual(result['decisions'][0]['decision'],'INSUFFICIENT_EVIDENCE')
            self.assertTrue(any(g['name']=='folding_ddg' and g['status']=='MISSING' for g in result['decisions'][0]['gates']))
