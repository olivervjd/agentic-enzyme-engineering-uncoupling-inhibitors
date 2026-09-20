import importlib.util
import json
from pathlib import Path
import unittest

from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.examples.prepare_epsps_chemical_states import enumerate_states

BASELINE = Path(__file__).resolve().parents[1] / 'herbicide_desensitization_agent/app/registry/epsps_evidence.json'


class EPSPSBaselineTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(BASELINE.read_text())

    def test_context_preserves_s3p_and_other_targets_metals(self):
        registry = TargetRegistry.default()
        epsps = registry.get('AT2G45300', 'glyphosate')
        self.assertEqual(epsps.required_context, ['shikimate-3-phosphate'])
        self.assertIn('Mg2+', registry.get('AT3G48560', 'chlorsulfuron').required_context)
        self.assertIn('divalent-metal-ions', registry.get('AT1G66200', 'glufosinate').required_context)

    def test_experimental_states_and_missing_neutral_are_explicit(self):
        structure = self.data['experimental_structures'][0]
        self.assertEqual(structure['pdb_id'], '7PXY')
        self.assertTrue(structure['same_target'])
        self.assertEqual(structure['state'], 'open_apo')
        self.assertFalse(structure['glyphosate_present'])
        neutral = next(c for c in self.data['controls'] if c['class'] == 'neutral')
        self.assertEqual(neutral['status'], 'UNAVAILABLE')
        self.assertEqual(self.data['conditions']['at2022_activity']['pH'], 7.4)
        self.assertEqual(self.data['conditions']['at2022_crystallization']['pH'], 7.5)

    def test_kinetics_are_never_direct_binding_or_invented_intervals(self):
        for measurement in self.data['measurements']:
            self.assertNotEqual(measurement['endpoint'], 'Kd')
            self.assertFalse(measurement['direct_binding_measurement'])
            self.assertIsNone(measurement['confidence_interval'])
            self.assertIn(measurement['conditions_id'], self.data['conditions'])
        tips_vmax = next(m for m in self.data['measurements'] if m['subject']=='ecoli_T97I_P101S' and m['endpoint']=='Vmax')
        self.assertEqual(tips_vmax['value'], 7)
        self.assertEqual(tips_vmax['reported_plus_minus'], .2)

    def test_canonical_mapping_is_not_homolog_numbering(self):
        mapping = {r['ecoli_residue']: r for r in self.data['mapping']['residues']}
        self.assertEqual(mapping['D313']['arabidopsis_canonical_residue'], 'D407')
        self.assertEqual(mapping['G96']['arabidopsis_canonical_residue'], 'G177')
        self.assertEqual(mapping['G96']['arabidopsis_mature_model_position'], 101)
        self.assertEqual([p['canonical'] for p in self.data['protected_catalytic_positions']], [99,407,435])

    @unittest.skipUnless(importlib.util.find_spec('dimorphite_dl') and importlib.util.find_spec('rdkit'), 'optional chemistry dependencies unavailable')
    def test_actual_enumeration_retains_uncertainty_and_stereochemistry(self):
        data = enumerate_states({'pep':'C=C(C(=O)O)OP(=O)(O)O', 's3p':'C1[C@H]([C@@H]([C@@H](C=C1C(=O)O)OP(=O)(O)O)O)O', 'glyphosate':'C(C(=O)O)NCP(=O)(O)O'})
        self.assertFalse(data['protonation_experimentally_validated'])
        for name, charge in [('pep',-3),('s3p',-3),('glyphosate',-2)]:
            ligand=data['ligands'][name]
            self.assertEqual(ligand['selected_campaign_state']['formal_charge'], charge)
            self.assertGreater(len(ligand['microstates']), 1)
            for state in ligand['microstates']:
                self.assertIsNone(state['population'])
                self.assertTrue(state['covalent_identity_and_stereochemistry_retained'])
        with self.assertRaises(ValueError):
            enumerate_states({'glyphosate':'CC'})
