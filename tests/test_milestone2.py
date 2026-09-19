import json
import tempfile
import unittest
from pathlib import Path

from herbicide_desensitization_agent.app.agents.pipeline_agents import InteractionFingerprintAgent
from herbicide_desensitization_agent.app.backends.bionemo import BioNeMoIRStructureBackend
from herbicide_desensitization_agent.app.backends.diffdock import DiffDockNIMBackend
from herbicide_desensitization_agent.app.chemistry.contacts import residue_contacts_from_pdb_and_sdf
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import Ligand, Pose, Provenance, StructureModel, TargetProtein
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.app.validators.structure_context_validator import validate_structure_context
from herbicide_desensitization_agent.app.visualization import MolstarArtifactRenderer


PROVENANCE = [Provenance("test://fixture", "unit-test", "synthetic")]


class FakeTransport:
    def __init__(self):
        self.calls = []

    def request(self, method, path, payload=None):
        self.calls.append((method, path, payload))
        if path == "/v1/health/ready":
            return {"status": "ready"}
        return {"poses": [{"confidence": 0.8, "contacts": [3, 7]}, {"confidence": 0.6, "contacts": [4]}]}


class Milestone2Tests(unittest.TestCase):
    def setUp(self):
        self.registry = TargetRegistry.default()

    def test_bionemo_processor_adapter_parses_ensemble(self):
        def request_factory(record_id, sequence, entry):
            return {"input_id": record_id, "sequence": sequence}

        def processor(rows):
            return [
                {
                    "__record_id": row["__record_id"],
                    "output_path": f"/tmp/{row['__record_id']}.cif",
                    "format": "cif",
                    "scores": json.dumps({"plddt": [80.0, 100.0], "ptm": 0.7}),
                }
                for row in rows
            ]

        entry = self.registry.get("AT2G45300", "glyphosate")
        backend = BioNeMoIRStructureBackend(
            processor,
            request_factory,
            ensemble_size=2,
            context_provider=lambda selected: list(selected.required_context),
        )
        models = backend.predict_ensemble(
            TargetProtein(entry.agi, entry.protein_name, "MALW", provenance=PROVENANCE), entry
        )
        self.assertEqual(len(models), 2)
        self.assertEqual(models[0].confidence, 0.9)
        self.assertEqual(models[0].format, "cif")

    def test_diffdock_nim_adapter_builds_pose_ensemble(self):
        with tempfile.TemporaryDirectory() as directory:
            protein = Path(directory) / "protein.pdb"
            protein.write_text("HEADER TEST\n", encoding="utf-8")
            structure = StructureModel(
                "model-1", "AT2G45300", "test", 0.8, [], PROVENANCE, str(protein), "pdb"
            )
            transport = FakeTransport()
            backend = DiffDockNIMBackend(transport, num_poses=2)
            poses = backend.dock([structure], Ligand("glyphosate", "herbicide", "C", provenance=PROVENANCE))
            self.assertEqual(len(poses), 2)
            self.assertEqual(poses[0].contacts, [3, 7])
            self.assertEqual(transport.calls[0][1], DiffDockNIMBackend.ENDPOINT)
            self.assertTrue(backend.healthcheck())

    def test_contact_fingerprint_uses_ensemble_sets(self):
        native = [Pose("n1", "m1", "native", 0.7, [2, 3], PROVENANCE)]
        herbicide = [Pose("h1", "m1", "herbicide", 0.7, [3, 5], PROVENANCE)]
        result = InteractionFingerprintAgent().compare("AT2G45300", native, herbicide, {7})
        self.assertEqual(result.shared_protected, [3])
        self.assertEqual(result.herbicide_selective_mutable, [5])
        self.assertEqual(result.native_ligand_critical_protected, [2])
        self.assertEqual(result.second_shell_candidates, [4, 6])

    def test_geometric_contact_extraction(self):
        pdb = (
            "ATOM      1  CA  ALA A   7       0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  ALA A   8      10.000  10.000  10.000  1.00 20.00           C\n"
        )
        sdf = "ligand\n  test\n\n  1  0  0  0  0  0            999 V2000\n   1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"
        self.assertEqual(residue_contacts_from_pdb_and_sdf(pdb, sdf), [7])

    def test_target_specific_context_validation(self):
        entry = self.registry.get("ATCG00020", "atrazine")
        isolated = StructureModel("m", entry.agi, "test", 0.8, ["protein-chain"], PROVENANCE)
        failures = validate_structure_context(isolated, entry)
        self.assertTrue(any("isolated soluble protein" in failure for failure in failures))

    def test_artifact_store_serializes_dataclass_lists(self):
        with tempfile.TemporaryDirectory() as directory:
            store = ArtifactStore(directory)
            path = store.write_manifest(
                "run-1", "poses.json", [Pose("p", "m", "ligand", 0.5, [1], PROVENANCE)]
            )
            self.assertEqual(json.loads(path.read_text())[0]["pose_id"], "p")

    def test_bionemo_adapter_reads_current_output_paths_envelope(self):
        row = {"output_paths": json.dumps({"cif": "/tmp/current.cif"})}
        self.assertEqual(BioNeMoIRStructureBackend._artifact_path(row), "/tmp/current.cif")

    def test_molstar_renderer_creates_standalone_viewer_bundle(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "prediction.cif"
            source.write_text("data_model\n_entry.id model\n", encoding="utf-8")
            model = StructureModel(
                "model-1", "AT2G45300", "boltz-2", 0.8, [], PROVENANCE,
                str(source), "cif", {"plddt": [0.7, 0.9], "ptm": 0.6},
            )
            renderer = MolstarArtifactRenderer(ArtifactStore(Path(directory) / "artifacts"))
            manifest = renderer.render("run-1", model, 1)
            run_dir = Path(directory) / "artifacts" / "run-1"
            viewer = (run_dir / manifest["viewer"]).read_text(encoding="utf-8")
            self.assertIn("molstar.Viewer.create", viewer)
            self.assertIn("Mean pLDDT: 0.800", viewer)
            self.assertTrue((run_dir / manifest["structure"]).is_file())
            self.assertTrue((run_dir / manifest["scores"]).is_file())


if __name__ == "__main__":
    unittest.main()
