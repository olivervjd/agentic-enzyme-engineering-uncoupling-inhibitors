"""Export fresh EPSPS calculations without inventing calibration or nominations."""
from __future__ import annotations
import argparse
from collections import defaultdict
import itertools
import json
import math
from pathlib import Path
import shutil
import statistics
import time

from .export_evidence_dashboard import (CONTACT_LEGEND, contact_records, conservation,
    replicated_geometry, experimental_geometry, table_export, read)
from .prepare_demo import geometry, inside
from .execution_table import execution_table
from ..app.backends.calibration import sha, write
from ..app.backends.structural import TMAlignStructuralMatcher
from ..app.orchestrator.execution import ExecutionJournal, utc_now

LIGAND_CONTEXT = {'glyphosate': 'herbicide', 'pep': 'native', 's3p': 'native'}
PROTECTED = {'arabidopsis_wt': {99,407,435}, 'ecoli_wt': {22,313,341}, 'ecoli_G96A': {22,313,341}}
BLOCKERS = [
    'Matched same-target resistant and catalytic-loss assay controls are incomplete; neutral controls are documented unavailable and are not an independent mandatory blocker.',
    'Seed ranges are descriptive; calibrated endpoint uncertainty and acceptance thresholds are unavailable.',
    'Published Km/Ki/IC50 kinetics cannot be substituted for Kd or model affinity scores.',
    'Affinity parser standardization can change ligand protonation; affinity geometry is excluded from charged-state contact analysis.',
    'Independent evidence requires matched chemical context and validated control recovery; external results alone do not establish this.',
    'No new experimental native-function or herbicide-resistance measurements were performed.',
    'No candidates reached a finalist stage, so molecular dynamics and alchemical free-energy calculations were not run.',
    'Genome-edit feasibility and legislative eligibility remain unevaluated: no candidate coding-sequence edit, crop jurisdiction or applicable framework has been supplied.',
]


def validate_predictions(manifest, predictions):
    expected = {(r['mode'],r['subject'],r['context'],r['seed']):r for r in manifest['records']}
    required = {(mode,subject,context,seed) for subject in manifest['subjects'] for seed in manifest['seeds']
                for mode,contexts in [('structure',('apo','s3p_only','herbicide','native')),
                                      ('affinity',('herbicide','native','native_s3p'))] for context in contexts}
    if len(manifest['seeds']) != 3 or set(manifest['subjects']) != {'arabidopsis_wt','ecoli_wt','ecoli_G96A'} or set(expected) != required:
        raise ValueError('Fresh campaign requires the complete preregistered 63-job panel')
    observed = [(r['mode'],r['subject'],r['context'],r['seed']) for r in predictions]
    if len(expected)!=len(manifest['records']) or len(set(observed))!=len(observed) or set(expected)!=set(observed):
        raise ValueError('Missing, extra or duplicate fresh prediction records')
    for row in predictions:
        original=expected[(row['mode'],row['subject'],row['context'],row['seed'])]
        for key in ('id','request','request_sha256','ligands'):
            if row[key]!=original[key]:
                raise ValueError('Prediction provenance differs from effective input manifest')
        if row['mode']=='affinity' and (isinstance(row.get('predicted_pIC50'),bool) or not isinstance(row.get('predicted_pIC50'),(int,float)) or not math.isfinite(row['predicted_pIC50'])):
            raise ValueError('Invalid affinity output')


def aggregate_contacts(raw, predictions, manifest, scores):
    """PEP and S3P independently share the same three native-complex seeds."""
    predictions=[p for p in predictions if p['mode']=='structure']
    structure_ids={p['id'] for p in predictions}
    relevant=[r for r in raw if r['pose_id'] in structure_ids and r['context']==LIGAND_CONTEXT[r['ligand']]]
    grouped=defaultdict(list)
    for r in relevant: grouped[(r['variant'],r['residue'],r['ligand'])].append(r)
    sites={(r['variant'],r['residue']) for r in relevant}
    sites.update((subject,position) for subject,positions in PROTECTED.items() if subject in manifest['subjects'] for position in positions)
    result=[]
    for subject,position in sorted(sites):
        counts={lig:len([p for p in predictions if p['subject']==subject and p['context']==ctx]) for lig,ctx in LIGAND_CONTEXT.items()}
        frequency={lig:len({r['pose_id'] for r in grouped[(subject,position,lig)]})/counts[lig] if counts[lig] else None for lig in LIGAND_CONTEXT}
        native=[frequency[lig] for lig in ('pep','s3p')]
        maximum=max(native) if all(v is not None for v in native) else None
        herb=frequency['glyphosate']; protected=position in PROTECTED.get(subject,set())
        classification='INSUFFICIENT_EVIDENCE'
        if protected: classification='CATALYTIC_OR_PROTECTED'
        elif maximum is not None and herb is not None:
            if maximum>=2/3: classification='SHARED_HERBICIDE_NATIVE_CONTACT' if herb>=2/3 else 'NATIVE_CRITICAL_CONTACT'
            elif herb>=2/3: classification='HERBICIDE_SELECTIVE_CONTACT'
            elif herb>0: classification='UNSTABLE_OR_METHOD_DEPENDENT_CONTACT'
        offset=manifest['subjects'][subject].get('offset',0)
        for ligand in LIGAND_CONTEXT:
            records=grouped[(subject,position,ligand)]
            nearest=min(records,key=lambda r:r['min_distance']) if records else None
            result.append({'variant':subject,'residue':position,'structure_residue':position-offset,'chain':'A','ligand':ligand,
                'amino_acid':manifest['subjects'][subject]['sequence'][position-offset-1],
                'context':LIGAND_CONTEXT[ligand],'model_mode':'structure','contact_frequency':frequency[ligand],
                'replicates':counts[ligand],'supporting_replicates':len({r['pose_id'] for r in records}),
                'min_distance':nearest['min_distance'] if nearest else None,'distance_unit':'angstrom',
                'interaction_types':['general proximity contact'] if records else [],
                'pose_ids':[r['pose_id'] for r in records],
                'protein_atoms':sorted({a['protein_atom'] for r in records for a in r['atom_pairs']}),
                'ligand_atoms':sorted({a['ligand_atom'] for r in records for a in r['atom_pairs']}),
                'classification':classification,'protected_status':'Protected catalytic homolog-mapped residue' if protected else 'Conservatively held: native contact' if maximum is not None and maximum>=2/3 else 'Unreviewed; calibration gate closed',
                'herbicide_contact_frequency':herb,'native_contact_frequencies':dict(zip(('pep','s3p'),native)),
                'maximum_native_contact_frequency':maximum,
                'herbicide_selectivity_difference':herb-maximum if herb is not None and maximum is not None else None,
                'herbicide_native_contact_ratio':herb/maximum if herb is not None and maximum else None,
                'conservation':scores.get((subject,position)), 'independent_method_count':1 if records else 0,
                'methods_supporting':['Boltz-2 2.2.1'] if records else [],
                'contact_uncertainty':'Observed seed fraction, no calibrated interval',
                'method_agreement':'Independent matched-context contact validation not established',
                'native_function_importance':'Catalytic homolog evidence' if protected else 'Proximity alone does not establish catalytic importance',
                'distance_from_catalytic_residues':None,'mutation_feasibility':'Not eligible; calibration incomplete'})
    return result


def binding_rows(predictions):
    affinity=[p for p in predictions if p['mode']=='affinity']
    result=[]
    for subject,ligand in sorted({(p['subject'],p['ligand']) for p in affinity}):
        group=sorted([p for p in affinity if (p['subject'],p['ligand'])==(subject,ligand)],key=lambda p:p['seed'])
        values=[p['predicted_pIC50'] for p in group]
        paired=[]
        if subject=='ecoli_G96A':
            controls={p['seed']:p for p in affinity if p['subject']=='ecoli_wt' and p['ligand']==ligand and p['context']==group[0]['context']}
            if set(controls)!={p['seed'] for p in group}: raise ValueError('G96A delta requires matched WT seeds')
            for p in group:
                wt=controls[p['seed']]
                if p['ligands']!=wt['ligands']: raise ValueError('Paired affinity chemistry differs')
                paired.append({'seed':p['seed'],'mutant_id':p['id'],'wt_id':wt['id'],'mutant_minus_wt':p['predicted_pIC50']-wt['predicted_pIC50']})
        deltas=[p['mutant_minus_wt'] for p in paired]
        result.append({'variant':subject,'ligand':ligand,'role':'herbicide' if ligand=='glyphosate' else 'native substrate',
            'evidence_state':'Insufficient evidence','affinity_metric':'predicted pIC50','unit':'dimensionless (-log10 IC50/M)',
            'estimate':statistics.mean(values),'observed_range':[min(values),max(values)],'interval':None,
            'interval_description':'Observed seed range is not a confidence interval','replicates':len(group),'independent_methods':1,
            'wt_relative_change':statistics.mean(deltas) if deltas else None,'paired_seed_deltas':paired,
            'delta_observed_range':[min(deltas),max(deltas)] if deltas else None,
            'delta_direction':'Negative model pIC50 change means weaker modeled inhibitor response; not a measured kinetic or thermodynamic change',
            'replicate_results':[{'id':p['id'],'seed':p['seed'],'value':p['predicted_pIC50'],'affinity_sha256':p['affinity_sha256']} for p in group],
            'effective_ligands':group[0]['ligands'],'geometry_scope':'Excluded: affinity parser state is separate from charged structure-mode geometry',
            'decision':'INSUFFICIENT_EVIDENCE','limitations':BLOCKERS})
    return result


def annotate_experimental_binding(bindings, ligand_geometry):
    """Resolved crystal ligand occupancy is scoped separately from model gates."""
    for subject,code,ligands,complex_ligands in [('ecoli_wt','1G6S',('glyphosate','s3p'),['glyphosate','s3p']),
                                               ('ecoli_G96A','1MI4',('s3p',),['s3p'])]:
        validated=[g for g in ligand_geometry if g['subject']==subject and g['reference']==code]
        if not validated: continue
        for ligand in ligands:
            if not all(ligand in g['ligand_atom_mapping'] for g in validated):
                raise ValueError('Experimental binding annotation requires validated deposited ligand identity')
            row=next((b for b in bindings if b['variant']==subject and b['ligand']==ligand),None)
            if row is None: continue
            row['experimental_evidence']={'evidence_ladder_state':'Experimentally confirmed binding',
                'observation':'Ligand resolved in a crystallographically determined protein complex',
                'scope':'Deposited E. coli homolog crystal complex only','subject':subject,'ligand':ligand,
                'pdb_id':code,'complex_ligands':complex_ligands,'method':'X-ray crystallography',
                'direct_Kd_measurement':False,'Kd':None,'ligand_protonation':None,'assay_pH':None,
                'equivalence_to_current_model_context':'Not established; crystal occupancy does not validate selected protonation or modeled affinity context',
                'transfer_to_Arabidopsis':'No same-target experimental ligand complex inferred',
                'provenance':{'source_type':'experimental','source':'https://www.rcsb.org/structure/'+code,
                              'artifact':{'path':'references/'+code+'.cif','sha256':validated[0]['reference_sha256']}},
                'limitations':['Crystal occupancy establishes this ligand-bound structural observation, not a binding constant or resistance phenotype',
                               'For S3P, the crystal co-ligand context differs from the modeled PEP+S3P native complex']}
    return bindings


def structural_metrics(root, manifest, predictions, raw):
    matcher=TMAlignStructuralMatcher(); result=[]
    structure=[p for p in predictions if p['mode']=='structure']
    for subject,context in sorted({(p['subject'],p['context']) for p in structure}):
        group=[p for p in structure if (p['subject'],p['context'])==(subject,context)]
        info=manifest['subjects'][subject]
        site=sorted({r['structure_residue'] for r in raw if r['pose_id'] in {p['id'] for p in group}} | {p-info.get('offset',0) for p in PROTECTED.get(subject,set())})
        for left,right in itertools.combinations(group,2):
            metrics=matcher.compare(inside(root,left['structure']),inside(root,right['structure']),target_sequence=info['sequence'],active_site_residues=site)
            result.append({'subject':subject,'context':context,'query':left['id'],'reference':right['id'],'kind':'repeatability','model_mode':'structure','active_site_definition':'Union of observed matched-context contacts plus mapped catalytic residues',**metrics})
        reference=None; reference_offset=0
        if subject=='arabidopsis_wt' and context=='apo': reference=root/'references/7PXY.cif'; reference_offset=-76
        elif (subject,context) in [('ecoli_wt','herbicide'),('ecoli_G96A','s3p_only')]: reference=inside(root,info['reference'])
        if reference:
            for p in group:
                metrics=matcher.compare(inside(root,p['structure']),reference,target_sequence=info['sequence'],reference_offset=reference_offset,active_site_residues=site)
                result.append({'subject':subject,'context':context,'query':p['id'],'reference':'7PXY' if reference_offset else info['pdb'],
                    'kind':'experimental_accuracy_same_target_apo' if reference_offset else 'experimental_accuracy_homolog',
                    'model_mode':'structure','reference_sha256':sha(reference),'reference_offset':reference_offset,
                    'interpretation':'State-matched structural comparison, not proof of function or resistance',**metrics})
    for mutant in [p for p in structure if p['subject']=='ecoli_G96A']:
        wt=next(p for p in structure if p['subject']=='ecoli_wt' and p['context']==mutant['context'] and p['seed']==mutant['seed'])
        site=sorted({r['structure_residue'] for r in raw if r['pose_id'] in {mutant['id'],wt['id']}} | PROTECTED['ecoli_wt'])
        metrics=matcher.compare(inside(root,mutant['structure']),inside(root,wt['structure']),target_sequence=manifest['subjects']['ecoli_wt']['sequence'],mutation='G96A',active_site_residues=site)
        result.append({'subject':'ecoli_G96A','context':mutant['context'],'query':mutant['id'],'reference':wt['id'],'kind':'matched_mutant_vs_wt','model_mode':'structure','seed':mutant['seed'],**metrics})
    return result


def baseline_availability(manifest, comparisons, baseline, root, docking_summaries=()):
    """Validate that baseline evidence exists, not that its values pass calibration."""
    from ..app.backends.calibration import validate_a3m
    from ..app.evidence_system.schema import digest
    root=Path(root); sources=baseline.get('sources',{}); entries={}
    def fact(name,available,source,artifact,scope,source_type='computational'):
        entries[name]={'validated':bool(available),'validation_scope':scope,
                       'acceptance_calibrated':False,'provenance':{'source_type':source_type,'source':source,'artifact':artifact}}
    measured=[c for c in comparisons if c.get('kind')=='experimental_accuracy_same_target_apo'
              and c.get('subject')=='arabidopsis_wt' and c.get('context')=='apo'
              and isinstance(c.get('query_tm_score'),(float,int)) and c.get('reference_sha256')]
    fact('experimental_structure_comparison',bool(measured),TMAlignStructuralMatcher.method,
         {'comparisons':measured,'sha256':digest(measured)},'Same-target apo comparison calculated with mapped coordinates; accuracy acceptance threshold remains uncalibrated')
    site_sources={k:sources[k] for k in ('ec2009','at2022') if k in sources}
    fact('published_binding_site',bool(site_sources) and bool(baseline.get('mapping',{}).get('residues')),
         'Primary EPSPS structures and canonical sequence mapping',{'sources':site_sources,'baseline_sha256':digest(baseline)},
         'Published homolog binding-site evidence mapped to the target; no same-target ligand-bound experimental pose claimed','experimental')
    protected=baseline.get('protected_catalytic_positions',[])
    sequence=manifest.get('subjects',{}).get('arabidopsis_wt',{}).get('sequence',''); offset=manifest.get('subjects',{}).get('arabidopsis_wt',{}).get('offset',76)
    catalytic=bool(protected) and all(0 <= p['canonical']-offset-1 < len(sequence) and sequence[p['canonical']-offset-1]==p['residue'] and p['source'] in sources for p in protected)
    fact('catalytic_residues',catalytic,'Primary homolog catalytic mutagenesis and verified target sequence mapping',
         {'positions':protected,'baseline_sha256':digest(baseline)},'Catalytic residue identities and mapping verified; homolog functional transfer remains labeled','experimental')
    msa_artifacts=[]; msa_errors=[]
    for subject,info in manifest.get('subjects',{}).items():
        try:
            path=inside(root,info['msa']['path'])
            if sha(path)!=info['msa']['sha256']: raise ValueError('MSA digest mismatch')
            validate_a3m(path.read_text(),info['sequence'])
            msa_artifacts.append({'subject':subject,'path':str(path),'sha256':sha(path)})
        except (KeyError,FileNotFoundError,ValueError) as exc:
            msa_errors.append({'subject':subject,'error':str(exc)})
    fact('documented_msa',bool(msa_artifacts) and not msa_errors,'Hashed ColabFold query-aligned MSA inputs',
         {'files':msa_artifacts,'errors':msa_errors},'Input hashes and complete query sequence validated; no new sequence-search claim')
    fact('biological_assembly','at2022' in sources,'Arabidopsis primary study: size exclusion and PDBePISA monomer assignment',
         {'source':sources.get('at2022'),'curated_baseline_sha256':digest(baseline)},
         'Published monomeric target supports modeled monomer; no mutant assembly stability conclusion','experimental')
    repeats=[c for c in comparisons if c.get('kind')=='repeatability' and c.get('subject')=='arabidopsis_wt']
    expected_pairs=len(manifest.get('seeds',[]))*(len(manifest.get('seeds',[]))-1)//2
    repeated=expected_pairs>=3 and all(len([c for c in repeats if c['context']==context])==expected_pairs for context in ('apo','s3p_only','herbicide','native'))
    fact('repeated_structure_prediction',repeated,'Boltz-2 independent seeds; recalculated coordinate comparison',
         {'comparisons':repeats,'sha256':digest(repeats),'seeds':manifest.get('seeds',[])},
         'Repeated structures and comparisons present; reproducibility acceptance thresholds remain uncalibrated')
    docking=[d for d in docking_summaries if d.get('independent_docking_seed_count',0)>=3 and d.get('provenance',{}).get('sha256')]
    fact('repeated_docking',bool(docking),'AutoDock Vina completed seed ensembles',
         {'jobs':[{k:d[k] for k in ('source_record','context','independent_docking_seed_count','reference_scope','provenance')} for d in docking]},
         'Repeated docking performed on explicitly recorded receptor contexts; no calibrated redocking pass or affinity agreement claimed')
    entries['neutral_mutation_controls']={'unavailable':True,'search_provenance':sources,
         'required_for_this_gate':False,'note':'Documented unavailable neutral control is permitted by the live registry; resistant and catalytic-loss assay control requirements remain'}
    return entries


def evaluate_readiness(manifest, bindings, contacts, comparisons, baseline, root, docking_summaries=()):
    """Run actual core gates with unavailable bounds; never invent a protocol."""
    from ..app.evidence_system import calibrate_thresholds, classify_residues, evaluate_candidate
    from ..app.evidence_system.calibration import validate_live_registry
    from ..app.evidence_system.schema import digest
    registry={'glyphosate':'herbicide','pep':'native_substrate','s3p':'native_substrate'}
    performed={p['context'] for p in manifest.get('records',[]) if p['mode']=='structure' and p.get('structure')}
    complexes=[]
    for context,names in [('herbicide',['glyphosate+S3P']),('native',['PEP+S3P','S3P+PEP']),('apo',['apo']),('s3p_only',['substrate_bound_control'])]:
        if context in performed: complexes.extend(names)
    spec={'target_id':'AT2G45300','herbicide_id':'glyphosate','native_ligand_ids':['pep','s3p'],
          'required_evaluations':registry,'assembly_required':False,'metrics':{},
          'evaluated_complexes':complexes,
          'baseline_evidence':baseline_availability(manifest,comparisons,baseline,root,docking_summaries),
          'prospective_rules':manifest.get('prospective_rules',{})}
    wt=[{'seed':seed,'method':'Boltz-2','binding_measurements':{b['ligand']:{'metric':'pIC50','unit':'dimensionless',
         'value':next(v['value'] for v in b['replicate_results'] if v['seed']==seed)} for b in bindings if b['variant']=='arabidopsis_wt' and b.get('affinity_metric')=='predicted pIC50'},
         'structural_comparisons':[c for c in comparisons if c['subject']=='arabidopsis_wt'],
         'provenance':{'source_type':'computational','source':'Boltz-2','artifact':str(root/'predictions.json')}} for seed in manifest['seeds']]
    try:
        core=calibrate_thresholds(spec,wt,baseline['controls'],registered_at=manifest['registered_at'],candidates_started_at=utc_now())
    except ValueError as exc:
        observed={b['ligand']:{'min':b['observed_range'][0],'max':b['observed_range'][1],
                  'range':b['observed_range'][1]-b['observed_range'][0],'n':b['replicates'],'unit':b['unit'],
                  'scope':'Descriptive affinity seed variability; not calibrated uncertainty'}
                  for b in bindings if b['variant']=='arabidopsis_wt' and b.get('observed_range')}
        core={'status':'INSUFFICIENT_EVIDENCE','specification':spec,'thresholds':{},'wt_variability':observed,
              'wt_replicates':wt,'input_controls':baseline['controls'],'limitations':[str(exc)],
              'registered_at':manifest['registered_at'],'binding':{},'contact':{}}
        core['digest']=digest(core)
    errors=validate_live_registry(core,target_id='AT2G45300',herbicide='glyphosate',native_ligands=['pep','s3p'],epsps=True)
    raw_errors=list(errors)
    if core.get('wt_variability'):
        errors=['Preregistered acceptance thresholds missing; descriptive WT seed variability is recorded' if e=='thresholds or WT variability missing' else e for e in errors]
    same_target_controls=[c for c in baseline.get('controls',[]) if c.get('scope')=='same_target']
    for kind in ('resistant','catalytic_loss'):
        if not any(c.get('class')==kind and c.get('status')!='UNAVAILABLE' for c in same_target_controls):
            errors.append('Matched same-target '+kind.replace('_','-')+' assay control missing')
    if any(b.get('interval') is None for b in bindings if b.get('affinity_metric')=='predicted pIC50'):
        errors.append('Calibrated model-affinity uncertainty intervals are unavailable')
    errors=list(dict.fromkeys(core.get('limitations',[])+errors))
    return {'core_calibration':core,'errors':errors,'eligible':not errors,'required_evaluations':registry,
            'raw_core_registry_errors':raw_errors,
            'prospective_threshold_status':'No complete preregistered acceptance/uncertainty protocol available; not inferred from observed WT results'}


def nomination_readiness(contacts, readiness):
    from ..app.evidence_system import classify_residues
    grouped=defaultdict(list)
    for r in contacts:
        if r['variant']=='arabidopsis_wt': grouped[r['residue']].append(r)
    rows=[]
    for position,items in sorted(grouped.items()):
        records={r['ligand']:{'frequency':r['contact_frequency'],'interval':None,
                  'method_frequencies':{'Boltz-2':r['contact_frequency']},'n':r['replicates']} for r in items}
        rows.append({'canonical_residue':position,'herbicide':records['glyphosate'],'native':{k:records[k] for k in ('pep','s3p')},
                     'conservation':items[0]['conservation'],'protected':position in PROTECTED['arabidopsis_wt'],
                     'catalytic':position in PROTECTED['arabidopsis_wt']})
    try:
        classified=classify_residues(rows,native_ligands=['pep','s3p'],thresholds=readiness['core_calibration'].get('contact',{}))
        errors=[]
    except ValueError as exc:
        classified=[]; errors=[str(exc)]
    eligible=[r for r in classified if r['nomination_eligible']] if readiness['eligible'] else []
    return {'input_residues':rows,'core_classification':classified,'eligible_residues':eligible,
            'errors':readiness['errors']+errors,'candidate_generation_enabled':bool(eligible) and readiness['eligible']}


def function_diagnostic(bindings, readiness, comparisons=()):
    from ..app.evidence_system import evaluate_candidate
    registry=readiness['required_evaluations']; records=[]
    for b in bindings:
        if b['variant']!='ecoli_G96A': continue
        wt=next(x for x in bindings if x['variant']=='ecoli_wt' and x['ligand']==b['ligand'])
        def affinity(row):
            return {'metric':'pIC50','unit':'dimensionless','value':row['estimate'],'interval':None,
                    'provenance':{'source_type':'computational','source':'Boltz-2','artifact':row['replicate_results']}}
        records.append({'evaluation_id':b['ligand'],'ligand':b['ligand'],'role':registry[b['ligand']],
            'wt_affinity':affinity(wt),'mutant_affinity':affinity(b),
            'wt_context':{'ligands':wt['effective_ligands']},'mutant_context':{'ligands':b['effective_ligands']},
            'assessment':{'variant':'ecoli_G96A','ligand':b['ligand'],'evaluation_id':b['ligand'],'methods':[]}})
    metrics={}
    for target,source,unit,aggregate in [('tm_score_forward','query_tm_score','dimensionless',min),
            ('tm_score_reverse','target_tm_score','dimensionless',min),('ca_lddt','alignment_lddt','dimensionless',min),
            ('coverage','alignment_coverage','dimensionless',min),('global_ca_rmsd','global_ca_rmsd_angstrom','angstrom',max),
            ('active_site_backbone_rmsd','active_site_rmsd_angstrom','angstrom',max)]:
        measured=[c for c in comparisons if c['kind']=='matched_mutant_vs_wt' and c.get(source) is not None]
        if measured:
            metrics[target]={'value':aggregate(c[source] for c in measured),'unit':unit,
                'provenance':{'source_type':'computational','source':TMAlignStructuralMatcher.method,'artifact':[{'query':c['query'],'reference':c['reference'],'value':c[source]} for c in measured]},
                'aggregation':'Conservative observed range endpoint across matched control contexts; not a calibrated interval'}
    result=evaluate_candidate({'variant':'ecoli_G96A','required_evaluations':registry,'bindings':records,'structural_metrics':metrics},readiness['core_calibration'])
    result.update(scope='Published homolog control diagnostic; not a new Arabidopsis nomination',model_binding_inputs=records,structural_metric_inputs=metrics)
    return result


def export(root, output, *, journal=None, stability=None, docking=None):
    started=time.monotonic(); root=Path(root); output=Path(output)
    if (output/'data/evidence_system.json').exists(): raise FileExistsError('Refusing to overwrite an existing evidence export')
    manifest=read(root/'effective_model_inputs.json'); predictions=read(root/'predictions.json')
    if isinstance(predictions,dict): predictions=predictions['records']
    validate_predictions(manifest,predictions)
    analysis=ExecutionJournal(output/'analysis_journal.json')
    contact_event=analysis.start('contact_mapping',inputs={'predictions_sha256':sha(root/'predictions.json'),'mode':'structure'},
                                 input_artifacts=[root/'predictions.json'],tools=['Biopython','numpy','MSA query identity'])
    from Bio.PDB import MMCIFParser
    parser=MMCIFParser(QUIET=True); models={}; raw=[]
    for p in predictions:
        for key,hashkey in [('structure','structure_sha256'),('request','request_sha256')]:
            if sha(inside(root,p[key]))!=p[hashkey]: raise ValueError('Changed prediction artifact: '+p['id'])
        if p['mode']=='affinity':
            apath=inside(root,p['structure']).parent/f"affinity_{p['id']}.json"
            if sha(apath)!=p['affinity_sha256']: raise ValueError('Changed affinity artifact')
        else:
            model=parser.get_structure(p['id'],inside(root,p['structure']))[0]; models[p['id']]=model
            raw.extend(contact_records(model,p,manifest['subjects'][p['subject']].get('offset',0)))
    structure=[p for p in predictions if p['mode']=='structure']
    contacts=aggregate_contacts(raw,structure,manifest,conservation(root,manifest))
    # Catalytic proximity uses actual heavy atoms in each matched charged structure.
    import numpy as np
    for row in contacts:
        distances=[]; offset=manifest['subjects'][row['variant']].get('offset',0)
        for p in structure:
            if p['subject']!=row['variant'] or p['context']!=row['context']: continue
            chain=models[p['id']]['A']; residue=chain[row['structure_residue']]
            others=[a.coord for position in PROTECTED.get(row['variant'],set()) if position!=row['residue']
                    for a in chain[position-offset] if a.element not in {'H','D'}]
            atoms=[a.coord for a in residue if a.element not in {'H','D'}]
            if atoms and others:
                distances.append(float(np.linalg.norm(np.array(atoms)[:,None,:]-np.array(others)[None,:,:],axis=-1).min()))
        row['distance_from_catalytic_residues']=min(distances) if distances else None
        row['catalytic_distance_definition']='Minimum heavy-atom distance to another mapped catalytic residue across same-context structure seeds; angstrom'
        row['catalytic_distance_observed_range']=[min(distances),max(distances)] if distances else None
    analysis.finish(contact_event,outputs={'contacts':contacts,'raw_contacts':raw},evidence_status='AVAILABLE')
    comparisons=analysis.execute('structure',lambda:structural_metrics(root,manifest,structure,raw)+replicated_geometry(models,structure,raw)+experimental_geometry(root,manifest,models,structure,raw),
        inputs={'structure_ids':[p['id'] for p in structure]},tools=['TM-align','Biopython Kabsch','mapped CA lDDT'])
    from .epsps_experimental_ligands import experimental_ligand_geometry
    ligand_geometry=analysis.execute('structure',lambda:experimental_ligand_geometry(root,manifest,structure),
        inputs={'operation':'Chemically validated experimental ligand RMSD after protein-only fit'},
        input_artifacts=[root/'references/1G6S.cif',root/'references/1MI4.cif'],
        tools=['CCD graph mapping','RDKit symmetry-aware CalcRMS','protein CA Kabsch'])
    for result in ligand_geometry:
        matches=[c for c in comparisons if c['kind']=='experimental_pocket_accuracy_homolog'
                 and all(c[k]==result[k] for k in ('query','subject','context','reference'))]
        if len(matches)!=1: raise ValueError('Experimental ligand geometry lacks one matching pocket comparison')
        matches[0].update(ligand_rmsd_angstrom=result['ligand_rmsd_angstrom'],ligand_atom_mapping=result['ligand_atom_mapping'],
                          ligand_geometry_audit=result)
        matches[0].pop('ligand_rmsd_missing_reason',None)
        matches[0]['limitations']=list(dict.fromkeys(matches[0]['limitations']+result['limitations']))
    bindings=annotate_experimental_binding(binding_rows(predictions),ligand_geometry)
    baseline_path=Path(__file__).parents[1]/'app/registry/epsps_evidence.json'; baseline=read(baseline_path)
    external={}
    external_artifacts={}; mutation_rows=[]; decisions=[]; docking_controls=[]; docking_table=None
    for name,path in [('stability',stability),('docking',docking)]:
        if path: external[name]={'source':str(path),'sha256':sha(path),'result':read(path),'interpretation':'Preserved external evidence; does not automatically satisfy calibration or alter nomination decisions'}
    from .fresh_external_evidence import docking_tables, stability_tables
    if docking:
        docking_table,external_artifacts['docking']=analysis.execute('docking',lambda:docking_tables(docking,bindings),
            inputs={'operation':'Aggregate completed docking output; no new Vina invocation'},input_artifacts=[docking],
            evidence_status='AVAILABLE',tools=['Vina result projection; top-pose seed frequencies kept separate'])
        docking_controls=docking_table['controls']; external['docking']['table_projection']=docking_table
    readiness=analysis.execute('calibration',lambda:evaluate_readiness(manifest,bindings,contacts,comparisons,baseline,root,docking_table['summaries'] if docking_table else ()),
        inputs={'registered_at':manifest['registered_at'],'prospective_rules':manifest.get('prospective_rules',{})},
        input_artifacts=[root/'effective_model_inputs.json',baseline_path],evidence_status=lambda r:'VALIDATED' if r['eligible'] else 'INSUFFICIENT_EVIDENCE',
        tools=['calibrate_thresholds','validate_live_registry'])
    nomination=analysis.execute('mutation_design',lambda:nomination_readiness(contacts,readiness),
        inputs={'calibration':readiness,'contact_rows':contacts},evidence_status='INSUFFICIENT_EVIDENCE',tools=['classify_residues'])
    if not nomination['candidate_generation_enabled']:
        analysis.skip('mutation_design','SKIPPED_DEPENDENCY','Candidate generation blocked by the executed calibration/contact eligibility gates',inputs={'reasons':nomination['errors']})
    diagnostic=analysis.execute('function_retention',lambda:function_diagnostic(bindings,readiness,comparisons),
        inputs={'bindings':bindings,'calibration':readiness},evidence_status='INSUFFICIENT_EVIDENCE',tools=['evaluate_candidate'])
    if docking_table:
        contacts.extend(docking_table['contacts']); raw.extend(docking_table['pose_contacts'])
    if stability:
        stability_table,external_artifacts['stability']=analysis.execute('function_retention',lambda:stability_tables(stability,readiness,contacts),
            inputs={'operation':'Evaluate DDGun-scored existing hypotheses with core gates; no new mutant structure prediction'},
            input_artifacts=[stability],evidence_status='INSUFFICIENT_EVIDENCE',tools=['evaluate_candidate','DDGun3D result projection'])
        mutation_rows=stability_table['mutations']; decisions=stability_table['decisions']; external['stability']['table_projection']=stability_table
    analysis.finalize()
    agents,execution=execution_table(read(journal) if journal else None,journal_path=journal)
    execution['analysis_journal']={'path':str(analysis.path),'sha256':sha(analysis.path),'run_id':analysis.data['run_id']}
    if journal is None:
        agents,execution=execution_table(analysis.data,journal_path=analysis.path)
    controls=[{'control':'Arabidopsis WT repeated structures and same-target apo accuracy','expected_outcome':'State-matched reproducibility and experimental comparison',
        'predicted_outcome':[c for c in comparisons if c['subject']=='arabidopsis_wt'], 'pass_fail':'INSUFFICIENT_EVIDENCE',
        'experimental_agreement':'7PXY compared to apo only; calculated metrics retained without unvalidated thresholds','implication':'Structural agreement alone does not calibrate binding or native function'}]
    controls += [{'control':'E. coli G96A '+b['ligand'],'expected_outcome':'Known inhibitor/substrate tradeoff in homolog','predicted_outcome':b['paired_seed_deltas'],
        'replicate_variability':b['delta_observed_range'],'experimental_agreement':'NON_EQUIVALENT_ENDPOINTS','pass_fail':'INSUFFICIENT_EVIDENCE',
        'implication':'Model pIC50 differences are diagnostic; published kinetics retain their endpoint and organism'} for b in bindings if b['variant']=='ecoli_G96A']
    controls += [{'control':name,'predicted_outcome':None,'pass_fail':'NOT_EVALUATED','implication':'Missing required matched control evidence'} for name in ['Same-target resistant assay control','Same-target catalytic-loss assay control','Calibrated affinity uncertainty']]
    controls.append({'control':'Same-target neutral assay control','predicted_outcome':None,'pass_fail':'DOCUMENTED_UNAVAILABLE',
        'implication':'Optional where unavailable after a documented search; not an independent calibration blocker','search_provenance':baseline['sources']})
    controls.extend(docking_controls)
    calibration={'status':'VALIDATED' if readiness['eligible'] else 'INSUFFICIENT_EVIDENCE','locked':False,'thresholds':{},'threshold_status':readiness['prospective_threshold_status'],'candidate_design_enabled':nomination['candidate_generation_enabled'],'reasons':readiness['errors'],'prospective_rules':manifest.get('prospective_rules',{}),'core_gate_results':readiness}
    contexts=[{'id':p['id'],'variant':p['subject'],'context':p['context'],'mode':p['mode'],'seed':p['seed'],'chemical_state':p['ligands'],
        'ph':manifest['ph'],'ph_status':'Reported assay pH used as ligand preparation assumption, not simulated protonation equilibrium',
        'ionic_strength':None,'protein_sequence':manifest['subjects'][p['subject']]['sequence'],
        'residue_offset':manifest['subjects'][p['subject']].get('offset',0),'msa':manifest['subjects'][p['subject']]['msa'],
        'provenance':{'request_sha256':p['request_sha256'],'structure_sha256':p['structure_sha256']}} for p in predictions]
    plan=[{'step':1,'action':'Measure matched same-target WT, resistant and catalytic-loss controls, adding neutral controls when available, with native-substrate kinetics and inhibitor response and independent replicate uncertainty.','readout':'Km/kcat and Ki or IC50 in defined conditions','gate':'Do not translate kinetic quantities into Kd'},
          {'step':2,'action':'Validate independent docking and actual encoded microstates against matched control complexes; test retained charge alternatives if conclusions depend on them.','readout':'Pose/contact and control recovery','gate':'No cross-context agreement claims'},
          {'step':3,'action':'Freeze validated endpoint-specific thresholds before nominating conservative substitutions.','readout':'Traceable control-calibrated screen','gate':'Native function, catalytic protection, stability and herbicide evidence all required'}]
    ligands=[]
    from rdkit import Chem
    from rdkit.Chem import Draw
    for name,smiles in manifest['ligands'].items():
        mol=Chem.MolFromSmiles(smiles); drawer=Draw.MolDraw2DSVG(340,190); drawer.DrawMolecule(mol); drawer.FinishDrawing()
        ligands.append({'id':name,'name':name,'smiles':smiles,'canonical_identifier':Chem.MolToInchiKey(mol),'svg':drawer.GetDrawingText(),
            'chemical_state':{'formal_charge':Chem.GetFormalCharge(mol),'status':'Declared structure-mode hypothesis; affinity effective states recorded separately'},'experimental_references':baseline['sources']})
    tables={'binding':bindings,'residue_interactions':contacts,'mutation_selection':mutation_rows,'model_calibration':controls,'agent_execution':agents}
    legends={'binding':'Only affinity-mode outputs. Seed mean/range and paired-seed mutant-WT differences are descriptive, never confidence intervals, Kd or experimental resistance.',
        'residue_interactions':CONTACT_LEGEND+' Only structure-mode models contribute. PEP and S3P each use the same native complexes with separate contact votes.',
        'mutation_selection':'Existing diagnostic and historical substitutions only; no new nominations. DDGun3D folding delta-delta-G uses a custom uncalibrated MSA protocol; positive means destabilizing. No new matched mutant complexes or native-function assays were performed.',
        'model_calibration':'Numerical structure diagnostics are not automatically pass/fail criteria; missing controls and endpoint uncertainty remain blocking.',
        'agent_execution':'Execution comes only from recorded journal events. Missing journal events remain not_started. Execution success is separate from evidence validity.'}
    payload={'schema_version':'1.0','status':'INSUFFICIENT_EVIDENCE','project':{'target_id':'AT2G45300','name':'EPSPS','organism':'Arabidopsis thaliana','herbicide':'Glyphosate'},
        'objective':'Minimal mutations that reduce herbicide interference while preserving native function','tables':tables,'decisions':decisions,
        'calibration':calibration,'limitations':list(dict.fromkeys(BLOCKERS+manifest['limitations'])), 'context_records':contexts,
        'structural_comparisons':comparisons,'contact_definitions':legends['residue_interactions'],'ligands':ligands,'validation_plan':plan,
        'table_legends':legends,'legends':{'contacts':legends['residue_interactions'],'binding':legends['binding'],'structure':'Charged structure-mode replicas only; apo7PXY compared only to apo. No arbitrary acceptance thresholds.'},
        'baseline':baseline,'external_evidence':external,'execution':execution,'nomination_readiness':nomination,'control_function_assessments':[diagnostic],
        'conclusion':{'directly_measured':'Published Arabidopsis WT kinetics/open apo7PXY and homolog mutant controls; no new assays',
            'computationally_predicted':f'{len(structure)} fresh structure-mode models and {len(predictions)-len(structure)} separate affinity jobs; coordinate metrics and contacts recalculated',
            'methods_agree_on':'Independent equivalent-context validation not established','nominations':'No new nominations; required calibration evidence remains incomplete',
            'uncertain':BLOCKERS,'next_experiment':plan[0]}, 'export_runtime_seconds':time.monotonic()-started}
    data=output/'data'; data.mkdir(parents=True,exist_ok=True)
    viewer={}
    for p in structure:
        subject=p['subject']; wt_subject='ecoli_wt' if subject=='ecoli_G96A' else subject
        ref=next(x for x in structure if x['subject']==wt_subject and x['context']==p['context'] and x['seed']==p['seed'])
        anchor=min((x for x in structure if x['subject']==wt_subject and x['context']==p['context']),key=lambda x:x['seed'])
        item={k:p[k] for k in ('id','subject','context','seed','ligand')}; offset=manifest['subjects'][subject].get('offset',0)
        item.update(geometry(inside(root,p['structure']),offset,inside(root,anchor['structure'])),sha256=p['structure_sha256'],wt_reference_id=ref['id'],alignment_anchor_id=anchor['id'],
            model_mode='structure',predicted_pIC50=None,contacts={k:[x+offset for x in v] for k,v in p['contacts'].items()},protocol_scope='Fresh charged-state structure-only prediction')
        ligand_names={ligand['id']:ligand['name'] for ligand in p['ligands']}
        for ligand in item['ligands']: ligand['name']=ligand_names.get(ligand['chain'],ligand['chain'])
        viewer[p['id']]=item
    sources=[{'id':key,'title':s.get('title',s.get('finding',key)),'url':s['url'],'evidence_level':s['kind'],'sha256':s.get('local_artifact_sha256')} for key,s in baseline['sources'].items()]
    from ..app.registry.loader import TargetRegistry
    results={'target':{'id':'AT2G45300','name':'EPSPS','organism':'Arabidopsis thaliana','herbicide':'Glyphosate'},
        'targets':[{'id':e.agi,'name':e.protein_name,'herbicide':e.herbicide} for e in TargetRegistry.default().entries],
        'run':{'id':root.name,'target_id':'AT2G45300','mode':'Fresh calculations','predictions':len(predictions),'structure_predictions':len(structure),'comparisons':len(comparisons),'passed':None,'status':'INSUFFICIENT_EVIDENCE','candidates_advanced':0},
        'sources':sources,'summary':[],'comparisons':comparisons,'checks':{},'blocking_reasons':BLOCKERS,'thresholds':{},'archive':{'rows':[]},'workflow':[],'binding':{'status':'INSUFFICIENT_EVIDENCE'},'legends':payload['legends']}
    for name,value in [('evidence_system',payload),('results',results),('structures',viewer),('pose_contacts',raw),('decisions',decisions),('chemical_contexts',contexts),('structural_comparisons',comparisons)]: write(data/(name+'.json'),value)
    for name,rows in tables.items(): table_export(data,name,rows,legends[name])
    # Strict campaign-relative allowlist; all63 outputs preserved, never credential/runtime directories.
    paths={Path('effective_model_inputs.json'),Path('calibration_inputs.json'),Path('predictions.json')}
    for p in predictions:
        paths.update([Path(p['structure']),Path(p['request'])])
        if p['mode']=='affinity': paths.add(Path(p['structure']).parent/f"affinity_{p['id']}.json")
        confidence=Path(p['structure']).parent/f"confidence_{p['id']}_model_0.json"
        if inside(root,confidence).is_file(): paths.add(confidence)
    paths.update(p.relative_to(root) for p in (root/'references').glob('*') if p.is_file())
    paths.update(Path(subject['msa']['path']) for subject in manifest['subjects'].values())
    for name in ('selected_chemical_states.json','chemical_state_enumeration.json','model_execution_journal.json','preparation_journal.json'):
        if (root/name).is_file(): paths.add(Path(name))
    for rel in sorted(paths):
        source=inside(root,rel); destination=output/'raw/fresh'/rel; destination.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(source,destination)
    for name,path in [('stability',stability),('docking',docking),('execution_journal',journal)]:
        if path: destination=output/'raw'/(name+'.json'); destination.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,destination)
    for name,artifacts in external_artifacts.items():
        base=Path(stability if name=='stability' else docking).resolve().parent
        for source in artifacts:
            relative=source.resolve().relative_to(base)
            destination=output/'raw'/name/relative; destination.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(source,destination)
    shutil.copyfile(baseline_path,output/'raw/epsps_evidence.json')
    (output/'limitations.md').write_text('# Limitations\n\n'+'\n'.join('- '+x for x in payload['limitations'])+'\n')
    (output/'experimental-validation-plan.md').write_text('# Experimental validation plan\n\n'+'\n\n'.join(f"{p['step']}. {p['action']}" for p in plan)+'\n')
    write(output/'workflow-manifest.json',{'schema_version':'1.0','status':'INSUFFICIENT_EVIDENCE','calibration':calibration,'execution':execution,'source_manifest_sha256':sha(root/'effective_model_inputs.json'),
        'artifacts':[{'path':str(p.relative_to(output)),'sha256':sha(p)} for p in sorted(output.rglob('*')) if p.is_file()]})
    return payload


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path); parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--journal',type=Path); parser.add_argument('--stability',type=Path); parser.add_argument('--docking',type=Path)
    args=parser.parse_args(); payload=export(args.root,args.output,journal=args.journal,stability=args.stability,docking=args.docking)
    print(json.dumps({'status':payload['status'],'binding_rows':len(payload['tables']['binding']),'contact_rows':len(payload['tables']['residue_interactions']),'nominations':0}))


if __name__=='__main__': main()
