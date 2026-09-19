from __future__ import annotations

import argparse
import json
from pathlib import Path

from herbicide_desensitization_agent.app.schemas.models import Provenance, StructureModel
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.app.visualization import MolstarArtifactRenderer


def render(structure: Path, scores: Path, output: Path) -> Path:
    score_data = json.loads(scores.read_text(encoding="utf-8"))
    model = StructureModel(
        model_id="h200-live-smoke",
        target_agi="validation-example",
        method="nvidia-bionemo-ir:boltz-2",
        confidence=MolstarArtifactRenderer._mean_plddt(score_data) or 0.0,
        context=["protein-chain"],
        provenance=[
            Provenance(
                source="NVIDIA BioNeMo Inference Runtime",
                method="boltz-2",
                evidence_type="predicted-structure",
            )
        ],
        artifact_path=str(structure),
        format=structure.suffix.lstrip("."),
        scores=score_data,
    )
    store = ArtifactStore(output)
    renderer = MolstarArtifactRenderer(store)
    artifact = renderer.render("h200-live-smoke", model, 1)
    manifest = store.write_manifest("h200-live-smoke", "visualizations.json", [artifact])
    return manifest.parent / artifact["viewer"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a BioNeMo structure as a standalone Mol* viewer")
    parser.add_argument("structure", type=Path, help="BioNeMo CIF/mmCIF/PDB output")
    parser.add_argument("scores", type=Path, help="BioNeMo scores JSON")
    parser.add_argument("--output", type=Path, default=Path("artifacts"))
    args = parser.parse_args()
    print(render(args.structure, args.scores, args.output))


if __name__ == "__main__":
    main()
