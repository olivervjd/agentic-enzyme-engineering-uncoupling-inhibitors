import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.examples.evaluate_epsps_campaign import read_campaign


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        sequence = "A" * 520
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        self.manifest = {"target_agi": "AT2G45300", "modeled_domain": [77, 520],
                         "sequence": sequence, "sequence_sha256": digest, "records": []}
        for condition in ("native", "herbicide"):
            for replicate in (1, 2):
                identity = f"WT-{condition}-{replicate}"
                self.manifest["records"].append({"record_id": identity, "mutation": "WT", "condition": condition,
                                                  "replicate": replicate, "sequence_sha256": digest})
                (self.directory / f"{identity}.cif").touch()

    def read(self):
        (self.directory / "campaign_manifest.json").write_text(json.dumps(self.manifest))
        return read_campaign(self.directory)

    def test_incomplete_replicates_are_rejected(self):
        self.manifest["records"].pop()
        with self.assertRaisesRegex(ValueError, "replicate"):
            self.read()

    def test_wrong_sequence_hash_is_rejected(self):
        self.manifest["records"][0]["sequence_sha256"] = "wrong"
        with self.assertRaisesRegex(ValueError, "sequence hash"):
            self.read()

    def test_missing_coordinate_artifact_is_rejected(self):
        next(self.directory.glob("*.cif")).unlink()
        with self.assertRaisesRegex(ValueError, "Missing coordinates"):
            self.read()

    def test_complete_manifest_is_accepted(self):
        _, groups = self.read()
        self.assertEqual(len(groups[("WT", "native")]), 2)
