from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any


class ArtifactStore:
    """Run-scoped, JSON-manifest artifact persistence with path containment."""

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)

    def run_directory(self, run_id: str) -> Path:
        safe = "".join(char for char in run_id if char.isalnum() or char in "-_")
        if not safe or safe != run_id:
            raise ValueError("run_id must contain only letters, numbers, hyphens, and underscores")
        directory = (self.root / safe).resolve()
        if self.root not in directory.parents:
            raise ValueError("run directory escapes artifact root")
        directory.mkdir(parents=True, exist_ok=True)
        return directory

    def write_manifest(self, run_id: str, name: str, value: Any) -> Path:
        if not name.endswith(".json") or Path(name).name != name:
            raise ValueError("manifest name must be a simple .json filename")
        path = self.run_directory(run_id) / name
        payload = self._jsonable(value)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        return path

    @classmethod
    def _jsonable(cls, value: Any) -> Any:
        if is_dataclass(value):
            return cls._jsonable(asdict(value))
        if isinstance(value, dict):
            return {str(key): cls._jsonable(item) for key, item in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [cls._jsonable(item) for item in value]
        return value
