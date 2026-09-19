from ..registry.loader import TargetRegistry
from ..schemas.models import TargetRegistryEntry


def validate_pairing(registry: TargetRegistry, agi: str, herbicide: str) -> TargetRegistryEntry:
    return registry.get(agi, herbicide)

