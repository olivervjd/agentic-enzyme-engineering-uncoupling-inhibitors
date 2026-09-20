"""Synthetic molecular fixtures validate preparation contracts, not biology."""
import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from herbicide_desensitization_agent.examples.epsps_vina import (
    boltz_atom_names, molecule_with_coordinates, ccd_atom_map, validate_context,
    combine_rigid_receptor, jobs_from_manifest, cluster_poses, pose_molecule, pose_contacts, contact_frequencies, validate_protein_coordinates,
)


def fixture(smiles):
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol,randomSeed=42)!=0:raise RuntimeError('Fixture embedding failed')
    AllChem.UFFOptimizeMolecule(mol)
    names=boltz_atom_names(smiles)
    atoms={name:{'element':mol.GetAtomWithIdx(i).GetSymbol().upper(),
                 'xyz':tuple(mol.GetConformer().GetAtomPosition(i))} for i,name in names.items()}
    return atoms,names


def pdbqt_atom(serial,name,xyz=(0.,0.,0.),atom_type='C'):
    x,y,z=xyz
    return f'ATOM  {serial:5d} {name:>4s} LIG Q   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 {atom_type}'


class VinaPreparationTests(unittest.TestCase):
    def test_protein_coordinates_allow_exact_terminal_oxygen_name_swap_only(self):
        a={('A',' 427 ','O'):(0.,0.,0.),('A',' 427 ','OXT'):(1.,0.,0.)}
        b={('A',' 427 ','O'):(1.,0.,0.),('A',' 427 ','OXT'):(0.,0.,0.)}
        self.assertEqual(len(validate_protein_coordinates(a,b)),2)
        b[('A',' 427 ','O')]=(3.,0.,0.)
        with self.assertRaisesRegex(ValueError,'moved'):validate_protein_coordinates(a,b)
        with self.assertRaisesRegex(ValueError,'omitted'):validate_protein_coordinates(a,{})

    def test_coordinate_dictionary_order_does_not_change_graph_mapping(self):
        smiles='O=C([O-])C[NH2+]CP(=O)([O-])[O-]'
        atoms,names=fixture(smiles)
        shuffled=dict(reversed(list(atoms.items())))
        mol=molecule_with_coordinates(smiles,shuffled,names)
        self.assertEqual(Chem.GetFormalCharge(mol),-2)
        self.assertEqual(Chem.MolToSmiles(mol),Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        for i,name in names.items():
            self.assertLess((mol.GetConformer().GetAtomPosition(i)-Chem.rdGeometry.Point3D(*atoms[name]['xyz'])).Length(),1e-10)

    def test_missing_or_extra_cosubstrate_atom_is_rejected(self):
        smiles='C=C(OP(=O)([O-])[O-])C(=O)[O-]'
        atoms,names=fixture(smiles)
        lost=dict(atoms);lost.pop(next(iter(lost)))
        with self.assertRaisesRegex(ValueError,'every heavy atom'):
            molecule_with_coordinates(smiles,lost,names)
        extra={**atoms,'EXTRA':{'element':'O','xyz':(0.,0.,0.)}}
        with self.assertRaises(ValueError):molecule_with_coordinates(smiles,extra,names)

    def test_bond_geometry_and_element_mismatch_fail_closed(self):
        smiles='CCO';atoms,names=fixture(smiles)
        bad={k:dict(v) for k,v in atoms.items()};bad[names[0]]['xyz']=(99.,99.,99.)
        with self.assertRaisesRegex(ValueError,'bond'):molecule_with_coordinates(smiles,bad,names)
        bad={k:dict(v) for k,v in atoms.items()};bad[names[0]]['element']='N'
        with self.assertRaisesRegex(ValueError,'element'):molecule_with_coordinates(smiles,bad,names)

    def test_s3p_stereochemistry_reflection_is_rejected(self):
        smiles='O=C([O-])C1=C[C@@H](OP(=O)([O-])[O-])[C@@H](O)[C@H](O)C1'
        atoms,names=fixture(smiles)
        molecule_with_coordinates(smiles,atoms,names)
        reflected={k:{**v,'xyz':(-v['xyz'][0],v['xyz'][1],v['xyz'][2])} for k,v in atoms.items()}
        with self.assertRaisesRegex(ValueError,'stereochemistry'):
            molecule_with_coordinates(smiles,reflected,names)

    def test_native_targets_keep_the_other_registered_substrate(self):
        record={'mode':'structure','context':'native','ligands':[{'name':'pep','smiles':'C'},{'name':'s3p','smiles':'O'}]}
        self.assertEqual(validate_context(record,'pep'),{'s3p'})
        self.assertEqual(validate_context(record,'s3p'),{'pep'})
        record['ligands'].pop()
        with self.assertRaisesRegex(ValueError,'required context'):
            validate_context(record,'pep')

    def test_affinity_standardization_and_apo_are_not_silently_reused(self):
        record={'mode':'affinity','context':'herbicide','ligands':[]}
        with self.assertRaisesRegex(ValueError,'never affinity'):
            validate_context(record,'glyphosate')
        record['mode']='structure';record['context']='apo'
        with self.assertRaisesRegex(ValueError,'apo docking'):
            validate_context(record,'glyphosate')

    def test_rigid_receptor_retains_every_cosubstrate_atom_without_torsion_directives(self):
        protein=pdbqt_atom(1,'CA')+'\n'+pdbqt_atom(2,'O',atom_type='OA')+'\n'
        co='ROOT\n'+pdbqt_atom(1,'P1',atom_type='P')+'\nBRANCH 1 2\n'+pdbqt_atom(2,'O1',atom_type='OA')+'\nTORSDOF 1\n'
        combined,counts=combine_rigid_receptor(protein,{'s3p':co})
        self.assertEqual(counts['rigid_cosubstrate_atom_counts'],{'s3p':2})
        self.assertEqual(len(combined.splitlines()),4)
        self.assertEqual([int(line[6:11]) for line in combined.splitlines()],[1,2,3,4])
        self.assertNotIn('ROOT',combined);self.assertNotIn('BRANCH',combined)
        self.assertIn('P1',combined)
        with self.assertRaises(ValueError):combine_rigid_receptor(protein,{'s3p':'ROOT\nENDROOT'})

    def test_manifest_expands_two_native_ligands_and_excludes_affinity(self):
        base={'id':'native','subject':'unit','mode':'structure','context':'native','structure':'native.cif'}
        manifest={'records':[base,{**base,'id':'aff','mode':'affinity'},{**base,'id':'apo','context':'apo'}]}
        self.assertEqual([(r['id'],lig) for r,lig in jobs_from_manifest(manifest)],[('native','pep'),('native','s3p')])
        with self.assertRaises(ValueError):jobs_from_manifest(manifest,{'not-present'})

    def test_pose_rmsd_is_in_receptor_frame_without_alignment(self):
        atoms,names=fixture('CCO')
        mol=molecule_with_coordinates('CCO',atoms,names)
        translated=Chem.Mol(mol)
        conf=translated.GetConformer()
        for i in range(translated.GetNumAtoms()):
            x,y,z=conf.GetAtomPosition(i);conf.SetAtomPosition(i,(x+4,y,z))
        assignments,clusters=cluster_poses([mol,Chem.Mol(mol),translated],threshold=2.)
        self.assertEqual(assignments,[0,0,1]);self.assertEqual(clusters,[[0,1],[2]])

    def test_pose_mapping_rejects_missing_heavy_atom(self):
        atoms,names=fixture('CCO');mol=molecule_with_coordinates('CCO',atoms,names)
        lines=[pdbqt_atom(i+1,name,atoms[name]['xyz']) for i,name in names.items()]
        mapped=pose_molecule(mol,'\n'.join(lines));self.assertEqual(mapped.GetNumAtoms(),3)
        with self.assertRaisesRegex(ValueError,'incomplete'):
            pose_molecule(mol,'\n'.join(lines[:-1]))

    def test_contact_votes_count_independent_seeds_and_keep_atom_pairs(self):
        atoms,names=fixture('CCO');mol=molecule_with_coordinates('CCO',atoms,names)
        near=tuple(mol.GetConformer().GetAtomPosition(0))
        protein=[{'chain':'A','structure_residue':1,'residue':77,'insertion_code':'','amino_acid':'ALA',
                  'atoms':{'CA':{'element':'C','xyz':(near[0]+3.,near[1],near[2])}}}]
        common={'record_id':'synthetic-fixture','subject':'unit','context':'native','ligand':'pep'}
        poses=[{**common,'seed':211,'pose_rank':2,'pose_id':'211:2'},
               {**common,'seed':211,'pose_rank':3,'pose_id':'211:3'},
               {**common,'seed':223,'pose_rank':1,'pose_id':'223:1'}]
        contacts=pose_contacts(mol,poses[0],protein)+pose_contacts(mol,poses[1],protein)
        frequencies=contact_frequencies(contacts,poses)
        self.assertEqual(len(frequencies),1)
        self.assertEqual(frequencies[0]['independent_seed_frequency'],.5)
        self.assertEqual(frequencies[0]['top_pose_seed_frequency'],0.)
        self.assertAlmostEqual(frequencies[0]['pose_fraction'],2/3)
        self.assertEqual(contacts[0]['residue'],77)
        self.assertEqual(contacts[0]['interaction_types'],['general proximity contact'])
        self.assertTrue(contacts[0]['atom_pairs'][0]['protein_atom']=='CA')
        self.assertTrue(all(pair['distance_angstrom']<=5. for r in contacts for pair in r['atom_pairs']))

    def test_ccd_mapping_uses_graph_not_deposited_atom_order(self):
        # Synthetic CCD-style atom names intentionally differ from RDKit order.
        atoms,names=fixture('CCO');renamed={names[0]:'X3',names[1]:'X1',names[2]:'X2'}
        posed={renamed[name]:value for name,value in atoms.items()}
        cif='''data_fixture
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG X2 O
LIG X3 C
LIG X1 C
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG X3 X1 SING
LIG X1 X2 SING
'''
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'fixture.cif';path.write_text(cif)
            mapping,count=ccd_atom_map(path,'LIG','CCO',posed)
        self.assertEqual(mapping,{0:'X3',1:'X1',2:'X2'})
        self.assertEqual(count,1)


if __name__=='__main__':unittest.main()
