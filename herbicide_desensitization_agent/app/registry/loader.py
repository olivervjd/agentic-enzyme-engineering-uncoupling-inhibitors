from __future__ import annotations

import json
from pathlib import Path

from ..schemas.models import Provenance, TargetRegistryEntry


class TargetRegistry:
    """Loads the immutable project registry (JSON is valid YAML)."""

    def __init__(self, entries: list[TargetRegistryEntry], version: str) -> None:
        self.entries = entries
        self.version = version
        self._by_pair = {(e.agi.upper(), e.herbicide.casefold()): e for e in entries}

    @classmethod
    def default(cls) -> "TargetRegistry":
        path = Path(__file__).with_name("arabidopsis_targets.yaml")
        raw = json.loads(path.read_text(encoding="utf-8"))
        provenance = Provenance(
            source=str(path), method="curated-project-registry", evidence_type="configuration"
        )
        entries = [TargetRegistryEntry(provenance=[provenance], **item) for item in raw["entries"]]
        return cls(entries, raw["registry_version"])

    def get(self, agi: str, herbicide: str) -> TargetRegistryEntry:
        try:
            return self._by_pair[(agi.upper(), herbicide.casefold())]
        except KeyError as exc:
            raise ValueError(f"Unsupported AGI–herbicide pairing: {agi} / {herbicide}") from exc

