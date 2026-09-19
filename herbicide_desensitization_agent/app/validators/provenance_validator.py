from dataclasses import fields, is_dataclass
from typing import Any

from ..schemas.models import Provenance


def validate_provenance(value: Any) -> None:
    """Recursively require non-empty provenance on scientific dataclasses."""
    if isinstance(value, Provenance):
        if not value.source or not value.method or not value.evidence_type:
            raise ValueError("Incomplete provenance record")
        return
    if is_dataclass(value):
        if hasattr(value, "provenance") and not getattr(value, "provenance"):
            raise ValueError(f"Missing provenance: {type(value).__name__}")
        for item in fields(value):
            validate_provenance(getattr(value, item.name))
    elif isinstance(value, dict):
        for item in value.values():
            validate_provenance(item)
    elif isinstance(value, (list, tuple, set)):
        for item in value:
            validate_provenance(item)

