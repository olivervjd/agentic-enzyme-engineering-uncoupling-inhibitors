from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from .interfaces import StructurePredictionBackend
from ..schemas.models import Provenance, StructureModel, TargetProtein, TargetRegistryEntry


class BioNeMoIRStructureBackend(StructurePredictionBackend):
    """Adapter for BioNeMo IR's public `build_processor` callable.

    Construct with `from_runtime` on a GPU host, or inject a compatible processor
    in tests. The processor contract is `list[row] -> list[row]`.
    """

    def __init__(
        self,
        processor: Callable[[list[dict[str, Any]]], list[dict[str, Any]]],
        request_factory: Callable[[str, str, TargetRegistryEntry], Any],
        model_source: str = "boltz-2",
        ensemble_size: int = 2,
        context_provider: Callable[[TargetRegistryEntry], list[str]] | None = None,
    ) -> None:
        if ensemble_size < 2:
            raise ValueError("Milestone 2 requires a structure ensemble of at least two models")
        self.processor = processor
        self.request_factory = request_factory
        self.model_source = model_source
        self.ensemble_size = ensemble_size
        self.context_provider = context_provider or (lambda entry: ["protein-chain"])

    @classmethod
    def from_runtime(
        cls,
        output_dir: str | Path,
        model_source: str = "boltz-2",
        ensemble_size: int = 2,
        num_sampling_steps: int = 200,
    ) -> "BioNeMoIRStructureBackend":
        try:
            from bionemo_ir.data.schemas import InputRequest, MSARecord, Polymer
            from bionemo_ir.pipeline.processor.engine_proc import EngineProcessorConfig, build_processor
            from bionemo_ir.pipeline.stages.configs import FeatureGeneratorStageConfig, WriterStageConfig
        except ImportError as exc:
            raise RuntimeError(
                "BioNeMo IR is optional; install the NVIDIA-supported `bionemo-ir` wheel on a compatible GPU host"
            ) from exc

        output_path = str(Path(output_dir).resolve())
        config = EngineProcessorConfig(
            model_source=model_source,
            runtime_args={"num_sampling_steps": num_sampling_steps}
            if model_source in {"boltz-1", "boltz-2", "openfold3"}
            else {},
            feature_generator_stage=FeatureGeneratorStageConfig(init_context={"random_seed": 42}),
            writer_stage=WriterStageConfig(output_path=output_path, format="cif"),
        )

        def request_factory(record_id: str, sequence: str, entry: TargetRegistryEntry):
            return InputRequest(
                input_id=record_id,
                polymers=[
                    Polymer(
                        polymer_type="protein",
                        chain_id=["A1"],
                        sequence=sequence,
                        msas=[MSARecord(content=f">{record_id}\n{sequence}\n")],
                    )
                ],
            )

        return cls(build_processor(config), request_factory, model_source, ensemble_size)

    def predict_ensemble(self, target: TargetProtein, entry: TargetRegistryEntry) -> list[StructureModel]:
        rows = []
        for index in range(1, self.ensemble_size + 1):
            record_id = f"{target.agi}-{self.model_source}-{index}"
            request = self.request_factory(record_id, target.sequence, entry)
            rows.append({"record": request, "__record_id": record_id, "random_seed": 41 + index})
        output_rows = list(self.processor(rows))
        if len(output_rows) != self.ensemble_size:
            raise RuntimeError("BioNeMo processor returned an incomplete structure ensemble")
        models = []
        for row in output_rows:
            if row.get("__inference_error__"):
                raise RuntimeError(f"BioNeMo inference failed: {row['__inference_error__']}")
            artifact = row.get("output_path")
            if not artifact and not row.get("output_raw"):
                raise RuntimeError("BioNeMo output contains no structure artifact")
            scores = row.get("scores", {})
            if isinstance(scores, str):
                scores = json.loads(scores)
            confidence = self._confidence(scores)
            models.append(
                StructureModel(
                    model_id=str(row.get("__record_id")),
                    target_agi=target.agi,
                    method=f"nvidia-bionemo-ir:{self.model_source}",
                    confidence=confidence,
                    context=self.context_provider(entry),
                    provenance=[
                        Provenance(
                            source="NVIDIA BioNeMo Inference Runtime",
                            method=self.model_source,
                            evidence_type="predicted-structure",
                        )
                    ],
                    artifact_path=str(artifact) if artifact else None,
                    format=str(row.get("format", "cif")),
                    scores=scores,
                )
            )
        return models

    @staticmethod
    def _confidence(scores: dict[str, Any]) -> float:
        values = scores.get("plddt") or []
        if values:
            confidence = sum(float(value) for value in values) / len(values)
            return round(confidence / 100.0 if confidence > 1.0 else confidence, 4)
        for key in ("confidence_score", "ptm"):
            if key in scores and scores[key] is not None:
                return round(float(scores[key]), 4)
        return 0.0
