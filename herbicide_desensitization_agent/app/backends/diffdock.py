from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from .http import JsonTransport
from .interfaces import DockingBackend
from ..chemistry.contacts import residue_contacts_from_pdb_and_sdf
from ..schemas.models import Ligand, Pose, Provenance, StructureModel


class DiffDockNIMBackend(DockingBackend):
    ENDPOINT = "/molecular-docking/diffdock/generate"

    def __init__(
        self,
        transport: JsonTransport,
        num_poses: int = 10,
        response_adapter: Callable[[Any], list[dict[str, Any]]] | None = None,
    ) -> None:
        if num_poses < 2:
            raise ValueError("Pose ensembles require at least two poses")
        self.transport = transport
        self.num_poses = num_poses
        self.response_adapter = response_adapter or self._default_response_adapter

    def dock(self, structures: list[StructureModel], ligand: Ligand) -> list[Pose]:
        poses: list[Pose] = []
        for structure in structures:
            if not structure.artifact_path:
                raise ValueError(f"Structure {structure.model_id} has no readable artifact path")
            protein = Path(structure.artifact_path).read_text(encoding="utf-8")
            payload = {
                "protein": protein,
                "ligand": ligand.structure,
                "num_poses": self.num_poses,
            }
            response = self.transport.request("POST", self.ENDPOINT, payload)
            records = self.response_adapter(response)
            if not records:
                raise RuntimeError("DiffDock NIM returned no poses")
            for rank, record in enumerate(records, 1):
                confidence = float(record.get("confidence", record.get("score", 0.0)))
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
                        confidence=max(0.0, min(1.0, confidence)),
                        contacts=[int(value) for value in contacts],
                        provenance=[
                            Provenance(
                                source="NVIDIA DiffDock NIM",
                                method="diffdock-generate",
                                evidence_type="predicted-pose",
                            )
                        ],
                        artifact_path=record.get("artifact_path"),
                        rank=rank,
                        metadata={key: value for key, value in record.items() if key not in {"contacts"}},
                    )
                )
        return poses

    def healthcheck(self) -> bool:
        response = self.transport.request("GET", "/v1/health/ready")
        return bool(response)

    @staticmethod
    def _default_response_adapter(response: Any) -> list[dict[str, Any]]:
        if isinstance(response, list):
            return response
        if isinstance(response, dict):
            for key in ("poses", "results", "data"):
                value = response.get(key)
                if isinstance(value, list):
                    return value
        raise RuntimeError(
            "Unrecognized DiffDock response; inspect the deployed NIM's /docs schema and inject a response_adapter"
        )
