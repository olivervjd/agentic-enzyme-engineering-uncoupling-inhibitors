"""Project external scientific outputs into tables without pooling endpoints."""
import math
from pathlib import Path
import re
import statistics

from .export_evidence_dashboard import read
from .prepare_demo import inside
from ..app.backends.calibration import sha


def _finite(value):
    return not isinstance(value,bool) and isinstance(value,(float,int)) and math.isfinite(value)


def load_docking(path):
    path=Path(path).resolve(); data=read(path); jobs=[]; failures=[]; artifacts=[path]
    if 'records' in data:
        for record in data['records']:
            if record.get('execution_status')!='COMPLETED':
                failures.append(record); continue
            if not all(re.fullmatch(r'[A-Za-z0-9_.-]+',str(record.get(k,''))) for k in ('id','ligand')):
                raise ValueError('Unsafe docking job identity')
            relative=Path(record['id']+'--'+record['ligand'])/'docking_results.json'
            declared=Path(record['result'])
            # Remote absolute paths may be relocated only by the exact known job layout.
            if tuple(declared.parts[-2:])!=tuple(relative.parts):
                raise ValueError('Docking result does not match its declared job identity')
            result_path=inside(path.parent,relative)
            if record.get('originating_result_sha256') and sha(result_path)!=record['originating_result_sha256']:
                raise ValueError('Relocated docking result changed from its originating hash')
            job=read(result_path)
            if job['preparation']['record_id']!=record['id'] or job['preparation']['docked_ligand']!=record['ligand']:
                raise ValueError('Docking result identity mismatch')
            jobs.append((result_path,job))
    else:
        jobs=[(path,data)]
    for result_path,job in jobs:
        if job.get('status')!='COMPLETED' or not job.get('rows'):
            raise ValueError('Completed docking result requires poses')
        if any(not _finite(r.get('vina_score_kcal_mol')) or not _finite(r.get('reference_pose_rmsd_angstrom')) for r in job['rows']):
            raise ValueError('Docking has invalid scores or RMSD')
        keys=[(r['seed'],r['pose_rank']) for r in job['rows']]
        if len(keys)!=len(set(keys)): raise ValueError('Duplicate docking seed/pose')
        artifacts.append(result_path)
        allowed=['preparation.json','docking_journal.json','protein.pdb','receptor.pdbqt','ligand.pdbqt','pose_contacts.json','contact_frequencies.json']
        allowed += [f'poses-{seed}.{extension}' for seed in {r['seed'] for r in job['rows']} for extension in ('sdf','pdbqt')]
        for name in allowed:
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',name): raise ValueError('Unsafe docking artifact name')
            candidate=inside(result_path.parent,name)
            if candidate.is_file(): artifacts.append(candidate)
    return jobs,failures,sorted(set(artifacts))


def docking_tables(path, bindings):
    jobs,failures,artifacts=load_docking(path)
    controls=[]; contacts=[]; raw=[]; summaries=[]
    for source,job in jobs:
        rows=job['rows']; first=rows[0]; audit=job['preparation']; seeds={r['seed'] for r in rows}
        top=[r for r in rows if r['pose_rank']==1]
        if {r['seed'] for r in top}!=seeds: raise ValueError('Each docking seed needs a top pose')
        scores=[r['vina_score_kcal_mol'] for r in top]; rmsd=[r['reference_pose_rmsd_angstrom'] for r in top]
        summary={'source_record':first['record_id'],'context':first['context'],'reference_scope':first['reference_scope'],
            'method':'AutoDock Vina','quantity':'empirical pose score, not thermodynamic delta-G, Ki, Kd or pIC50',
            'unit':'kcal/mol','top_pose_score_mean':statistics.mean(scores),'top_pose_score_observed_range':[min(scores),max(scores)],
            'top_pose_reference_rmsd_observed_range_angstrom':[min(rmsd),max(rmsd)],
            'all_pose_reference_rmsd_observed_range_angstrom':[min(r['reference_pose_rmsd_angstrom'] for r in rows),max(r['reference_pose_rmsd_angstrom'] for r in rows)],
            'independent_docking_seed_count':len(seeds),'receptor_structure_count':1,'returned_pose_count':len(rows),
            'pose_clusters':job['pose_clusters'],'cluster_cutoff_angstrom':job['clustering_threshold_angstrom'],
            'cluster_cutoff_scope':'Pose clustering only; not a redocking acceptance threshold',
            'replicate_results':rows,'uncertainty_interval':None,
            'chemical_context_equivalent_to_boltz_affinity':False,
            'non_equivalence_reason':'Docked charged ligand and rigid co-substrate preparation differ from the affinity-standardized Boltz binder; do not claim independent affinity agreement',
            'method_assumptions':job['limitations'],'preparation':audit,
            'provenance':{'path':str(source),'sha256':sha(source)}}
        summaries.append(summary)
        binding=next((b for b in bindings if b['variant']==first['subject'] and b['ligand']==first['ligand']),None)
        if binding is None:
            binding={'variant':first['subject'],'ligand':first['ligand'],'role':'herbicide' if first['ligand']=='glyphosate' else 'native substrate',
                     'evidence_state':'Insufficient evidence','pose_hypothesis':'Docked poses available; binding support uncalibrated','affinity_metric':'No equivalent affinity evaluated','estimate':None,'interval':None,
                     'replicates':0,'independent_methods':0,'decision':'INSUFFICIENT_EVIDENCE'}
            bindings.append(binding)
        binding.setdefault('docking_evidence',[]).append(summary)
        binding.setdefault('pose_clusters',[]).append({'source_record':first['record_id'],'clusters':job['pose_clusters']})
        binding['method_agreement']='Independent pose diagnostics available; no equivalent-context calibrated affinity agreement'
        controls.append({'control':f"Vina {first['subject']} {first['ligand']} on {first['record_id']}",
            'expected_outcome':'Reproduce reference geometry under an independently validated redocking protocol',
            'predicted_outcome':summary,'replicate_variability':summary['top_pose_reference_rmsd_observed_range_angstrom'],
            'experimental_agreement':'Experimental pose comparison' if first['reference_scope']=='experimental_redocking' else 'Comparison to predicted reference pose',
            'pass_fail':'INSUFFICIENT_EVIDENCE','implication':'Geometry calculated; no preregistered redocking acceptance threshold. The 2 angstrom cutoff groups poses only.'})
        prefix=first['record_id']+'--'+first['ligand']+':'
        for row in job.get('pose_contacts',[]):
            raw.append({**row,'pose_id':prefix+row['pose_id'],'method':'AutoDock Vina','model_mode':'docking',
                        'provenance':summary['provenance']})
        for row in job.get('residue_contact_frequencies',[]):
            contacts.append({**row,'method':'AutoDock Vina','methods_supporting':['AutoDock Vina'],'independent_method_count':1,
                'model_mode':'docking','replicates':row['replicate_count'],'supporting_replicates':len(row['seeds_supporting_top_pose']),
                'pose_ids':[prefix+x for x in row['poses_supporting']], 'classification':'INSUFFICIENT_EVIDENCE',
                'protected_status':'Separate docking diagnostic; nomination eligibility not assessed',
                'mutation_feasibility':'Calibration incomplete','contact_uncertainty':row['uncertainty'],
                'method_agreement':'Separate receptor/ligand preparation; never pooled with Boltz contact votes',
                'frequency_scope':'Top-ranked pose per docking seed on this single receptor only',
                'provenance':summary['provenance']})
    for failure in failures:
        controls.append({'control':'Vina '+str(failure.get('id'))+' '+str(failure.get('ligand')),
            'pass_fail':'NOT_EVALUATED','execution_status':failure.get('execution_status'),
            'predicted_outcome':None,'implication':failure.get('reason','Docking computation not completed')})
    return {'summaries':summaries,'controls':controls,'contacts':contacts,'pose_contacts':raw,'failures':failures},artifacts


def stability_tables(path, readiness, contacts):
    from ..app.evidence_system import evaluate_candidate
    path=Path(path).resolve(); source=read(path); rows=[]; decisions=[]
    if source.get('method')!='DDGun3D' or source.get('unit')!='kcal/mol': raise ValueError('Unrecognized stability endpoint')
    seen=set()
    for mutation in source['mutations']:
        name=mutation['mutation']
        if name in seen or not _finite(mutation['folding_ddg_kcal_mol']): raise ValueError('Invalid or duplicate stability mutation')
        seen.add(name)
        provenance={'source_type':'computational','source':'DDGun3D '+source['protocol'],'artifact':{'path':str(path),'sha256':sha(path),'mutation':name}}
        result=evaluate_candidate({'variant':name,'required_evaluations':readiness['required_evaluations'],'bindings':[],
            'structural_metrics':{'folding_ddg':{'value':mutation['folding_ddg_kcal_mol'],'unit':'kcal/mol','provenance':provenance}}},readiness['core_calibration'])
        result.update(mutation=name,scope=mutation['scope'],new_nomination=False)
        decisions.append(result)
        residue=mutation['canonical_position']; site=next((r for r in contacts if r['variant']=='arabidopsis_wt' and r['residue']==residue),{})
        rows.append({'mutation':name,'residue':residue,'variant':'arabidopsis_wt','scope':mutation['scope'],
            'classification':result['decision'],'recommendation':result['decision'],'new_nomination':False,
            'selection_origin':'Homolog-mapped diagnostic position, not an experimentally verified Arabidopsis phenotype' if mutation['scope']=='known_control' else 'Historical exploratory hypothesis; not newly nominated',
            'ddg':mutation['folding_ddg_kcal_mol'],'ddg_unit':'kcal/mol','ddg_quantity':'Estimated folding delta-delta-G; positive means destabilizing',
            'native_ddgun_unfolding_ddg':mutation['ddgun3d_unfolding_ddg_kcal_mol'],'uncertainty_interval':mutation['uncertainty_interval'],
            'calibration_status':mutation['calibration_status'],'stability_protocol':source['protocol'],'stability_replicates':mutation['replicates'],
            'structure_retention':None,'binding_evidence':'No new matched Arabidopsis mutant complexes or affinity outputs were generated for this substitution',
            'herbicide_contacts_disrupted':None,'native_contacts_affected':None,'conservation':site.get('conservation'),
            'protected_status':site.get('protected_status'),'rationale':{'why_selected':'Previously defined diagnostic/historical panel',
                'expected_effect':'Not established by stability alone','native_interactions_to_retain':['PEP','S3P'],
                'falsification':'Measured catalytic loss or no inhibitor weakening invalidates a resistance hypothesis'},
            'limitations':source['limitations'],'gates':result['gates'],'provenance':provenance})
    # Only named local sibling products, with original hashes, enter the portable bundle.
    artifacts=[path]
    allowed={'7PXY_A_normalized_1_444.pdb','colabfold_query_aligned.fasta','wt.fasta','wt.colabfold_profile.hssp',
             '7PXY_A.dssp','ddgun3d_native.tsv','ddgun3d_features.json'}
    for item in source.get('output_artifacts',[]):
        name=Path(item['path']).name
        if name not in allowed: raise ValueError('Stability artifact is outside the scientific output allowlist')
        candidate=inside(path.parent,name)
        if not candidate.is_file() or sha(candidate)!=item['sha256']: raise ValueError('Missing or changed stability output artifact')
        artifacts.append(candidate)
    return {'mutations':rows,'decisions':decisions,'source':{'path':str(path),'sha256':sha(path)},'limitations':source['limitations']},sorted(set(artifacts))
