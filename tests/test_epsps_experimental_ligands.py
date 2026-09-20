import unittest
import numpy as np
from rdkit import Chem
from herbicide_desensitization_agent.examples.epsps_experimental_ligands import kabsch_transform,transformed_ligand_rmsd


def molecule(coords):
    mol=Chem.MolFromSmiles('CCO');conf=Chem.Conformer(3)
    for i,xyz in enumerate(coords):conf.SetAtomPosition(i,list(map(float,xyz)))
    mol.AddConformer(conf);return mol


class ExperimentalLigandTests(unittest.TestCase):
    def test_global_protein_transform_recovers_pose_without_ligand_fit(self):
        fixed=np.array([[0.,0.,0.],[3.,0.,0.],[0.,2.,0.],[0.,0.,4.]])
        q=np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]])
        shift=np.array([10.,20.,-7.]);moving=fixed@q+shift
        r,t,rmsd=kabsch_transform(moving,fixed)
        self.assertLess(rmsd,1e-10)
        pose=np.array([[1.,1.,1.],[2.,1.,1.],[2.,2.,1.]])
        self.assertLess(transformed_ligand_rmsd(molecule(pose@q+shift),molecule(pose),r,t),1e-6)
        # A rigid ligand translation remains four angstroms after protein fitting.
        displaced=(pose+np.array([4.,0.,0.]))@q+shift
        self.assertAlmostEqual(transformed_ligand_rmsd(molecule(displaced),molecule(pose),r,t),4.,places=6)

    def test_kabsch_rejects_invalid_matching(self):
        with self.assertRaises(ValueError):kabsch_transform([[0,0,0]],[[0,0,0]])
        with self.assertRaises(ValueError):kabsch_transform(np.zeros((3,3)),np.zeros((4,3)))
        with self.assertRaises(ValueError):kabsch_transform(np.full((3,3),np.nan),np.zeros((3,3)))

if __name__=='__main__':unittest.main()
