from __future__ import annotations

from pathlib import Path
from math import isfinite
from typing import Any, Callable

from .http import JsonTransport
from .interfaces import DockingBackend
from ..chemistry.contacts import residue_contacts_from_pdb_and_sdf
from ..schemas.models import Ligand, Pose, Provenance, StructureModel


class DiffDockNIMBackend(DockingBackend):
    ENDPOINT = "/molecular-docking/diffdock/generate"
    HOSTED_ENDPOINT = "/v1/biology/mit/diffdock"

    def __init__(
        self,
        transport: JsonTransport,
        num_poses: int = 10,
        response_adapter: Callable[[Any], list[dict[str, Any]]] | None = None,
        endpoint: str = ENDPOINT,
    ) -> None:
        if not 2 <= num_poses <= 100:
            raise ValueError("Pose ensembles require between 2 and 100 poses")
        self.transport = transport
        self.num_poses = num_poses
        self.response_adapter = response_adapter or self._default_response_adapter
        if endpoint not in {self.ENDPOINT, self.HOSTED_ENDPOINT}:
            raise ValueError("Unsupported DiffDock endpoint")
        self.endpoint = endpoint

    def dock(self, structures: list[StructureModel], ligand: Ligand) -> list[Pose]:
        poses: list[Pose] = []
        for structure in structures:
            if not structure.artifact_path:
                raise ValueError(f"Structure {structure.model_id} has no readable artifact path")
            raw_protein = Path(structure.artifact_path).read_text(encoding="utf-8")
            protein = "\n".join(line for line in raw_protein.splitlines() if line.startswith("ATOM  "))
            if not protein:
                raise ValueError("DiffDock requires protein PDB ATOM records; convert mmCIF before docking")
            formats = {"smiles": "txt", "txt": "txt", "sdf": "sdf", "mol2": "mol2"}
            if ligand.structure_format.lower() not in formats:
                raise ValueError("DiffDock ligand format must be SMILES, SDF, or MOL2")
            payload = {
                "protein": protein,
                "ligand": ligand.structure,
                "ligand_file_type": formats[ligand.structure_format.lower()],
                "num_poses": self.num_poses,
                "time_divisions": 20,
                "steps": 18,
                "save_trajectory": False,
            }
            response = self.transport.request("POST", self.endpoint, payload)
            records = self.response_adapter(response)
            if len(records) < 2:
                raise RuntimeError("DiffDock NIM returned fewer than two poses")
            for rank, record in enumerate(records, 1):
                if "confidence" not in record and "score" not in record:
                    raise ValueError("DiffDock pose is missing its confidence score")
                confidence = float(record.get("confidence", record.get("score")))
                if not isfinite(confidence):
                    raise ValueError("DiffDock returned non-finite pose confidence")
                contacts = record.get("contacts")
                if contacts is None:
                    sdf = record.get("ligand_sdf") or record.get("sdf") or record.get("pose")
                    contacts = residue_contacts_from_pdb_and_sdf(
                        protein, sdf, protein_chain=structure.scores.get("target_chain"),
                        residue_offset=structure.scores.get("structure_numbering_offset", 0),
                    ) if isinstance(sdf, str) else []
                poses.append(
                    Pose(
                        pose_id=f"{structure.model_id}-{ligand.role}-{rank}",
                        model_id=structure.model_id,
                        ligand_name=ligand.name,
                        confidence=confidence,
                        contacts=[int(value) for value in contacts],
                        provenance=[
                            Provenance(
                                source="NVIDIA DiffDock NIM",
                                method="diffdock-generate",
                                evidence_type="predicted-pose",
                                notes="Raw pose confidence, not probability or affinity. Receptor excludes HETATM context.",
                            )
                        ],
                        artifact_path=record.get("artifact_path"),
                        rank=rank,
                        metadata={**{key: value for key, value in record.items() if key != "contacts"},
                                  "confidence_kind": "raw-diffdock-pose-score", "receptor_context": "protein-only"},
                    )
                )
        return poses

    def healthcheck(self) -> bool:
        if self.endpoint == self.HOSTED_ENDPOINT:
            raise NotImplementedError("Hosted DiffDock has no local readiness endpoint; use an inference probe")
        response = self.transport.request("GET", "/v1/health/ready")
        return response is True or (isinstance(response, dict) and response.get("status") == "ready")

    @staticmethod
    def _default_response_adapter(response: Any) -> list[dict[str, Any]]:
        if isinstance(response, list):
            return response
        if isinstance(response, dict):
            if "status" in response and response["status"] != "success":
                raise RuntimeError("DiffDock prediction did not report success")
            if "ligand_positions" in response:
                positions = response["ligand_positions"]
                confidence = response.get("position_confidence")
                if not isinstance(positions, list) or not isinstance(confidence, list) or len(positions) != len(confidence):
                    raise ValueError("DiffDock pose and confidence arrays must have matching lengths")
                if not all(isinstance(sdf, str) and sdf.strip() for sdf in positions):
                    raise ValueError("DiffDock returned an empty or invalid pose")
                return [{"ligand_sdf": sdf, "confidence": score} for sdf, score in zip(positions, confidence)]
            for key in ("poses", "results", "data"):
                value = response.get(key)
                if isinstance(value, list):
                    return value
        raise RuntimeError(
            "Unrecognized DiffDock response; inspect the deployed NIM's /docs schema and inject a response_adapter"
        )
