from ..registry.loader import TargetRegistry
from ..schemas.models import WorkflowRequest
from .provenance_validator import validate_provenance
from .target_registry_validator import validate_pairing

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def validate_request(request: WorkflowRequest, registry: TargetRegistry):
    entry = validate_pairing(registry, request.target.agi, request.herbicide.name)
    sequence = request.target.sequence
    if not sequence or set(sequence) - AA20:
        raise ValueError("Protein sequence must contain only the 20 canonical amino acids")
    if request.target.residue_numbering_start != 1:
        raise ValueError("Use full sequence numbering starting at 1 and explicit structure offsets")
    if request.target.name != entry.protein_name:
        raise ValueError("Target protein name does not match the fixed registry")
    native_names = {ligand.name.casefold() for ligand in request.native_ligands}
    missing_ligands = {name.casefold() for name in entry.native_ligands} - native_names
    if missing_ligands:
        raise ValueError(f"Missing native ligands: {sorted(missing_ligands)}")
    context = {partner.name.casefold() for partner in request.functional_context}
    missing_context = {name.casefold() for name in entry.required_context} - context
    if missing_context:
        raise ValueError(f"Missing required functional context: {sorted(missing_context)}")
    validate_provenance(request)
    return entry
