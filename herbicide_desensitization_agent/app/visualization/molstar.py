from __future__ import annotations

import base64
import html
import json
from pathlib import Path
from typing import Any

from ..schemas.models import StructureModel
from ..storage.artifact_store import ArtifactStore


class MolstarArtifactRenderer:
    """Create a standalone Mol* view from a BioNeMo structure artifact."""

    MOLSTAR_VERSION = "4.18.0"

    def __init__(self, store: ArtifactStore) -> None:
        self.store = store

    def render(self, run_id: str, model: StructureModel, index: int) -> dict[str, str]:
        if not model.artifact_path:
            raise ValueError("structure model has no artifact path to visualize")
        source = Path(model.artifact_path)
        structure_format = (model.format or source.suffix.lstrip(".") or "cif").lower()
        if structure_format not in {"cif", "mmcif", "pdb"}:
            raise ValueError(f"Mol* visualization does not support {structure_format!r}")

        stem = f"structure-{index}"
        structure_name = f"{stem}.{structure_format}"
        score_name = f"{stem}-scores.json"
        viewer_name = f"{stem}-viewer.html"
        structure_path = self.store.copy_artifact(run_id, structure_name, source)
        score_path = self.store.write_manifest(run_id, score_name, model.scores)
        html = self._html(
            title=f"{model.target_agi} — {model.model_id}",
            structure_text=structure_path.read_text(encoding="utf-8"),
            structure_format="mmcif" if structure_format in {"cif", "mmcif"} else "pdb",
            scores=model.scores,
        )
        viewer_path = self.store.write_text(run_id, viewer_name, html)
        return {
            "model_id": model.model_id,
            "structure": structure_path.name,
            "scores": score_path.name,
            "viewer": viewer_path.name,
        }

    def _html(self, title: str, structure_text: str, structure_format: str, scores: dict[str, Any]) -> str:
        encoded = base64.b64encode(structure_text.encode("utf-8")).decode("ascii")
        mean_plddt = self._mean_plddt(scores)
        safe_title = json.dumps(title)
        safe_format = json.dumps(structure_format)
        confidence = "n/a" if mean_plddt is None else f"{mean_plddt:.3f}"
        ptm = self._formatted_score(scores.get("ptm"))
        iptm = self._formatted_score(scores.get("iptm"))
        version = self.MOLSTAR_VERSION
        return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@{version}/build/viewer/molstar.css">
  <style>
    html, body {{ margin: 0; height: 100%; font-family: system-ui, sans-serif; background: #111827; color: white; }}
    body {{ display: grid; grid-template-rows: minmax(240px, 1fr) auto auto; }}
    #viewer {{ position: relative; min-height: 240px; }}
    #summary {{ display: flex; flex-wrap: wrap; gap: 12px;
      align-items: center; padding: 12px 16px; box-sizing: border-box; background: #111827; font-size: 13px; }}
    #legend {{ margin: 0; padding: 0 16px 12px; font-size: 12px; line-height: 1.5; overflow-wrap: anywhere; }}
    #error {{ color: #fca5a5; }}
  </style>
</head>
<body>
  <div id="viewer"></div>
  <div id="summary"><strong id="title"></strong><span>Mean pLDDT: {confidence}</span><span>pTM: {ptm}</span><span>ipTM: {iptm}</span><span id="error"></span></div>
  <p id="legend">Figure legend: Predicted coordinates for the named model. Mean pLDDT describes local
  model confidence, pTM describes predicted global accuracy, and ipTM describes interface confidence;
  scores are displayed on a 0-1 scale, with higher values indicating greater confidence. These are
  model-confidence scores, not measured structural similarity, binding affinity, or retained function.
  Colors follow the selected Mol* representation. Missing scores are shown as n/a.</p>
  <script src="https://cdn.jsdelivr.net/npm/molstar@{version}/build/viewer/molstar.js"></script>
  <script>
    const title = {safe_title};
    document.getElementById('title').textContent = title;
    const structureData = atob('{encoded}');
    molstar.Viewer.create('viewer', {{ layoutIsExpanded: false, layoutShowControls: true }}).then(async viewer => {{
      await viewer.loadStructureFromData(structureData, {safe_format});
    }}).catch(error => {{ document.getElementById('error').textContent = error.message; }});
  </script>
</body>
</html>
"""

    @staticmethod
    def _mean_plddt(scores: dict[str, Any]) -> float | None:
        values = scores.get("plddt") or []
        if values:
            mean = sum(float(value) for value in values) / len(values)
        elif scores.get("mean_plddt") is not None:
            mean = float(scores["mean_plddt"])
        else:
            return None
        return mean / 100.0 if mean > 1.0 else mean

    @staticmethod
    def _formatted_score(value: Any) -> str:
        return "n/a" if value is None else f"{float(value):.3f}"
