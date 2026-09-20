"""Ligand geometry against homolog crystal poses, after a protein-only fit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from .epsps_vina import (boltz_atom_names, ccd_atom_map, molecule_with_coordinates,
                         heavy_atoms, select_ligand, digest, save)
from .prepare_demo import inside


def kabsch_transform(moving, fixed):
    import numpy as np
    moving=np.asarray(moving,dtype=float);fixed=np.asarray(fixed,dtype=float)
    if moving.shape!=fixed.shape or moving.ndim!=2 or moving.shape[1]!=3 or len(moving)<3:
        raise ValueError('Matched protein coordinates with at least three CA atoms required')
    if not np.isfinite(moving).all() or not np.isfinite(fixed).all():raise ValueError('Nonfinite protein coordinates')
    mc=moving.mean(axis=0);fc=fixed.mean(axis=0)
    u,_,vt=np.linalg.svd((moving-mc).T@(fixed-fc))
    if np.linalg.det(u@vt)<0:u[:,-1]*=-1
    rotation=u@vt;translation=fc-mc@rotation
    rmsd=float(np.sqrt(np.mean(np.sum((moving@rotation+translation-fixed)**2,axis=1))))
    return rotation,translation,rmsd


def transformed_ligand_rmsd(moving,fixed,rotation,translation):
    """Symmetry-aware RMSD in fitted receptor frame; never fit ligand atoms."""
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    result=Chem.Mol(moving);conf=result.GetConformer()
    for i in range(result.GetNumAtoms()):
        xyz=np.asarray(tuple(conf.GetAtomPosition(i)))@rotation+translation
        conf.SetAtomPosition(i,xyz.tolist())
    return float(rdMolAlign.CalcRMS(result,fixed,maxMatches=10000))


def experimental_ligand_geometry(root,manifest,predictions):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    root=Path(root);parser=MMCIFParser(QUIET=True);rows=[]
    references={('ecoli_wt','herbicide'):('1G6S',{'glyphosate':'GPJ','s3p':'S3P'}),
                ('ecoli_G96A','s3p_only'):('1MI4',{'s3p':'S3P'})}
    for record in predictions:
        key=(record['subject'],record['context'])
        if record.get('mode')!='structure' or key not in references:continue
        code,components=references[key];reference_path=root/'references'/f'{code}.cif'
        predicted_path=inside(root,record['structure'])
        predicted=parser.get_structure(record['id'],predicted_path)[0];reference=parser.get_structure(code,reference_path)[0]
        protein=[r for r in reference['A'] if r.id[0]==' ' and 'CA' in r]
        query=[r for r in predicted['A'] if r.id[0]==' ' and 'CA' in r]
        registered=manifest['subjects'][record['subject']]
        if [r.id[1] for r in protein]!=registered['observed_author_residue_numbers'] or [r.id[1] for r in query]!=list(range(1,len(protein)+1)):
            raise ValueError('Protein author-to-model numbering mapping differs from registered observed sequence')
        if ''.join(seq1(r.resname) for r in protein)!=registered['sequence'] or ''.join(seq1(r.resname) for r in query)!=registered['sequence']:
            raise ValueError('Protein sequence differs from registered matched reference')
        rotation,translation,ca_rmsd=kabsch_transform([r['CA'].coord for r in query],[r['CA'].coord for r in protein])
        values={};audits={}
        if {s['name'] for s in record['ligands']}!=set(components):raise ValueError('Experimental ligand context mismatch')
        for spec in record['ligands']:
            if spec.get('standardized_by_affinity_parser'):raise ValueError('Affinity-standardized coordinates are not supported')
            smiles=spec.get('model_smiles') or spec['smiles'];name=spec['name']
            predicted_atoms=heavy_atoms(select_ligand(predicted,spec));predicted_mapping=boltz_atom_names(smiles)
            experimental_atoms=heavy_atoms(select_ligand(reference,{'name':name,'component':components[name],'chain':'A'}))
            reference_mapping,mapping_count=ccd_atom_map(reference_path,components[name],smiles,experimental_atoms)
            pmol=molecule_with_coordinates(smiles,predicted_atoms,predicted_mapping)
            rmol=molecule_with_coordinates(smiles,experimental_atoms,reference_mapping)
            values[name]=transformed_ligand_rmsd(pmol,rmol,rotation,translation)
            audits[name]={'predicted_atom_map':predicted_mapping,'experimental_atom_map':reference_mapping,
                'equivalent_graph_mappings':mapping_count,'heavy_atom_count':pmol.GetNumAtoms(),'modeled_smiles':smiles,
                'mapping':'CCD element/bond-order graph isomorphism plus coordinate stereochemistry; Boltz canonical-rank names'}
        rows.append({'subject':record['subject'],'context':record['context'],'query':record['id'],'reference':code,
            'kind':'experimental_ligand_geometry_homolog','model_mode':'structure','seed':record['seed'],
            'ligand_rmsd_angstrom':values,'ligand_atom_mapping':audits,'global_ca_fit_rmsd_angstrom':ca_rmsd,
            'protein_ca_atom_count':len(query),'protein_mapping':'Registered observed sequence and author numbering; chain A',
            'protein_transform':{'rotation':rotation.tolist(),'translation':translation.tolist(),'convention':'row-vector coordinates @ rotation + translation'},
            'method':'Global protein CA Kabsch; RDKit symmetry-aware heavy-atom CalcRMS, no ligand fit',
            'predicted_sha256':digest(predicted_path),'reference_sha256':digest(reference_path),
            'limitations':['Homolog geometry control, not Arabidopsis ligand-pose validation',
                'Experimental ligand protonation is not established by coordinates; selected modeled charge is a mapping hypothesis, not confirmed chemical equivalence',
                'Geometric pose agreement does not establish binding affinity or resistance; acceptance threshold remains uncalibrated']})
    return rows


def main():
    from ..app.orchestrator.execution import ExecutionJournal
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign-dir',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.campaign_dir.resolve();manifest_path=root/'effective_model_inputs.json'
    manifest=json.loads(manifest_path.read_text());predictions=manifest['records']
    journal=ExecutionJournal(args.output.with_suffix('.execution.json'))
    try:
        rows=journal.execute('structure',lambda:experimental_ligand_geometry(root,manifest,predictions),
            inputs={'operation':'Experimental homolog ligand RMSD after protein-only Kabsch'},input_artifacts=[manifest_path],evidence_status='AVAILABLE')
        save(args.output,rows)
    finally:journal.finalize()
    print(json.dumps({'models':len(rows),'ligand_metrics':sum(len(r['ligand_rmsd_angstrom']) for r in rows)}))

if __name__=='__main__':main()
