import json
from pathlib import Path
import tempfile
import unittest

from herbicide_desensitization_agent.app.backends.calibration import sha
from herbicide_desensitization_agent.examples.package_evidence_dashboard import figure_scope_rows, package
from herbicide_desensitization_agent.examples.review_evidence_bundle import fingerprint


class PackageScopeTests(unittest.TestCase):
    def test_figures_do_not_pool_docking_or_omit_same_target_apo(self):
        boltz={'variant':'arabidopsis_wt','residue':99,'model_mode':'structure','contact_frequency':1.}
        docking={**boltz,'model_mode':'docking','method':'AutoDock Vina','contact_frequency':0.}
        affinity={'variant':'arabidopsis_wt','affinity_metric':'predicted pIC50','replicate_results':[{'value':1.}]}
        pose_only={'variant':'arabidopsis_wt','affinity_metric':'No equivalent affinity evaluated'}
        crystal={'kind':'experimental_accuracy_same_target_apo','reference':'7PXY'}
        binding,contacts,experimental=figure_scope_rows({'tables':{'binding':[affinity,pose_only],'residue_interactions':[boltz,docking]},'structural_comparisons':[crystal]})
        self.assertEqual(binding,[affinity]); self.assertEqual(contacts,[boltz]); self.assertEqual(experimental,[crystal])

    def test_second_embed_preserves_figure_bytes_and_rejects_wrong_review_content(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'bundle'; dist=Path(folder)/'dist'; dist.mkdir(); (root/'data').mkdir(parents=True); (root/'figures').mkdir()
            figure=root/'figures/fixed.svg'; figure.write_text('stable byte artifact')
            (dist/'index.html').write_text('<html><head></head><body></body></html>')
            evidence={'status':'INSUFFICIENT_EVIDENCE'}
            for name,value in [('evidence_system',evidence),('results',{}),('structures',{})]:
                (root/'data'/f'{name}.json').write_text(json.dumps(value))
            (root/'workflow-manifest.json').write_text('{}')
            source=root/'data/evidence_system.json'; original=source.read_bytes()
            review={'source_bundle_sha256':sha(source),'source_content_sha256':fingerprint(evidence),'source_snapshot_current':True,'status':'COMPLETED'}
            review_path=root/'data/evidence_model_review.json'; review_path.write_text(json.dumps(review))
            package(root,dist,render_figures=False)
            self.assertEqual(figure.read_text(),'stable byte artifact'); self.assertEqual(source.read_bytes(),original)
            self.assertIn('window.__EVIDENCE_DATA__',(root/'report.html').read_text())
            review['source_content_sha256']='wrong-content-digest'; review_path.write_text(json.dumps(review))
            with self.assertRaises(ValueError): package(root,dist,render_figures=False)
