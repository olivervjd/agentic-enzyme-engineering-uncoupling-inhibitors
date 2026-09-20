import json
import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.app.backends.calibration import sha, write
from herbicide_desensitization_agent.examples.epsps_fresh_campaign import prepare, SEEDS


class FreshCampaignTest(unittest.TestCase):
    def test_structure_and_affinity_are_distinct_protocols_with_complete_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            donor = base / 'donor'
            (donor / 'msa').mkdir(parents=True)
            (donor / 'references').mkdir()
            sequence = 'ACDEFGHIKLMNPQRSTVWY'
            msa = donor / 'msa/query.a3m'
            msa.write_text('>query\n' + sequence + '\n>homolog\nCCDEFGHIKLMNPQRSTVWY\n')
            write(donor / 'calibration_inputs.json', {'sequence': sequence,
                'subjects': {name: {'sequence': sequence, 'msa': {'path': 'msa/query.a3m', 'sha256': sha(msa)}}
                             for name in ('arabidopsis_wt', 'ecoli_wt', 'ecoli_G96A')}})
            reference = base / '7PXY.cif'
            reference.write_text('data_reference\n')
            root = base / 'fresh'
            data = prepare(root, donor, reference, remote_root='/model/campaign')
            self.assertEqual(len(data['records']), 63)
            self.assertFalse(data['prediction_reuse'])
            self.assertEqual(len({r['id'] for r in data['records']}), 63)
            for row in data['records']:
                request = json.loads((root / row['request']).read_text())
                self.assertEqual(sha(root / row['request']), row['request_sha256'])
                self.assertIn(row['seed'], SEEDS)
                self.assertEqual(request['sequences'][0]['protein']['msa'], '/model/campaign/msa/query.a3m')
                if row['mode'] == 'structure':
                    self.assertNotIn('properties', request)
                else:
                    binder = 'C' if row['context'] == 'native_s3p' else 'B'
                    self.assertEqual(request['properties'], [{'affinity': {'binder': binder}}])
                if row['context'] in {'native', 'native_s3p', 'herbicide'}:
                    self.assertEqual(len(request['sequences']), 3)
                    self.assertIn('[O-]', request['sequences'][-1]['ligand']['smiles'])
            with self.assertRaises(FileExistsError):
                prepare(root, donor, reference)


if __name__ == '__main__':
    unittest.main()
