import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.app.backends.diffdock import DiffDockNIMBackend
from herbicide_desensitization_agent.app.schemas.models import Ligand, StructureModel


PDB = "ATOM      1  CA  ALA A   7       0.000   0.000   0.000  1.00 20.00           C\n"
SDF = "ligand\n  test\n\n  1  0  0  0  0  0            999 V2000\n   1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"


class Transport:
    def __init__(self, response):
        self.response = response

    def request(self, method, path, payload=None):
        self.path, self.payload = path, payload
        return self.response


class DiffDockLiveContractTests(unittest.TestCase):
    def test_current_response_preserves_negative_confidence_and_contacts(self):
        transport = Transport({"status": "success", "ligand_positions": [SDF, SDF],
                               "position_confidence": [1.2, -2.7]})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "protein.pdb"
            path.write_text(PDB)
            structure = StructureModel("wt", "AT2G45300", "test", .8, [], [], str(path), "pdb",
                                       {"target_chain": "A", "structure_numbering_offset": 76})
            backend = DiffDockNIMBackend(transport, 2, endpoint=DiffDockNIMBackend.HOSTED_ENDPOINT)
            poses = backend.dock([structure], Ligand("glyphosate", "herbicide", "C"))
        self.assertEqual([pose.confidence for pose in poses], [1.2, -2.7])
        self.assertEqual(poses[0].contacts, [83])
        self.assertEqual(poses[0].metadata["confidence_kind"], "raw-diffdock-pose-score")
        self.assertEqual(transport.path, DiffDockNIMBackend.HOSTED_ENDPOINT)

    def test_rejects_mismatched_arrays_and_failed_status(self):
        with self.assertRaises(ValueError):
            DiffDockNIMBackend._default_response_adapter({"ligand_positions": [SDF], "position_confidence": []})
        with self.assertRaises(RuntimeError):
            DiffDockNIMBackend._default_response_adapter({"status": "failed", "ligand_positions": []})

    def test_healthcheck_does_not_accept_truthy_error_object(self):
        self.assertFalse(DiffDockNIMBackend(Transport({"status": "not ready"})).healthcheck())
        self.assertTrue(DiffDockNIMBackend(Transport(True)).healthcheck())

    def test_rejects_header_only_receptor(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.pdb"
            path.write_text("HEADER no atoms\n")
            structure = StructureModel("wt", "AT2G45300", "test", .8, [], [], str(path))
            with self.assertRaisesRegex(ValueError, "ATOM"):
                DiffDockNIMBackend(Transport({})).dock([structure], Ligand("x", "herbicide", "C"))

    def test_rejects_missing_or_nonfinite_confidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "protein.pdb"
            path.write_text(PDB)
            structure = StructureModel("wt", "AT2G45300", "test", .8, [], [], str(path))
            for records in ([{"contacts": []}, {"contacts": []}],
                            [{"confidence": float("nan"), "contacts": []}] * 2):
                with self.subTest(records=records), self.assertRaises(ValueError):
                    DiffDockNIMBackend(Transport(records)).dock([structure], Ligand("x", "herbicide", "C"))

    def test_hosted_readiness_is_not_faked(self):
        backend = DiffDockNIMBackend(Transport(True), endpoint=DiffDockNIMBackend.HOSTED_ENDPOINT)
        with self.assertRaises(NotImplementedError):
            backend.healthcheck()


if __name__ == "__main__":
    unittest.main()
