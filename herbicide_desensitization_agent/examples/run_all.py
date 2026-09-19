from __future__ import annotations

import json

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import (
    FunctionalPartner,
    Ligand,
    Provenance,
    TargetProtein,
    WorkflowRequest,
)


PLACEHOLDER_SEQUENCE = "MALWMRLLPLLALLALWGPGPGAGGAGGAGGAG"
PLACEHOLDER = [
    Provenance(
        source="example://placeholder",
        method="synthetic-fixture",
        evidence_type="synthetic",
        notes="Not an Arabidopsis protein or chemical structure",
    )
]


def make_request(entry) -> WorkflowRequest:
    return WorkflowRequest(
        target=TargetProtein(entry.agi, entry.protein_name, PLACEHOLDER_SEQUENCE, provenance=PLACEHOLDER),
        herbicide=Ligand(entry.herbicide, "herbicide", "C", provenance=PLACEHOLDER),
        native_ligands=[Ligand(name, "native", "C", provenance=PLACEHOLDER) for name in entry.native_ligands],
        functional_context=[
            FunctionalPartner(name, "required-context", provenance=PLACEHOLDER) for name in entry.required_context
        ],
        protected_residues={1, 2},
    )


def main() -> None:
    registry = TargetRegistry.default()
    backend = MockScientificBackend()
    workflow = WorkflowOrchestrator(registry, backend, backend, backend, backend, backend)
    summary = []
    for entry in registry.entries:
        result = workflow.run(make_request(entry))
        summary.append(
            {
                "agi": entry.agi,
                "herbicide": entry.herbicide,
                "status": result.packets[0].status,
                "candidate": result.packets[0].candidate.mutation,
                "synthetic": True,
            }
        )
    print(json.dumps({"runs": summary}, indent=2))


if __name__ == "__main__":
    main()

