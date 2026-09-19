from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from herbicide_desensitization_agent.app.backends.mocks import MockScientificBackend
from herbicide_desensitization_agent.app.backends.precomputed import PrecomputedComplexBackend
from herbicide_desensitization_agent.app.orchestrator.workflow import WorkflowOrchestrator
from herbicide_desensitization_agent.app.registry.loader import TargetRegistry
from herbicide_desensitization_agent.app.schemas.models import (
    FunctionalPartner, Ligand, Pose, Provenance, StructureModel, TargetProtein, WorkflowRequest,
)
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore


UNIPROT_PROVENANCE = [Provenance(
    "https://rest.uniprot.org/uniprotkb/P05466", "UniProt reviewed record", "curated-sequence-and-sites",
)]
PUBCHEM_PROVENANCE = [Provenance(
    "https://pubchem.ncbi.nlm.nih.gov", "PubChem PUG REST", "curated-chemical-record",
)]
BIONEMO_PROVENANCE = [Provenance(
    "NVIDIA BioNeMo Inference Runtime", "boltz-2 protein-ligand complex, 100 steps", "predicted-structure",
    "Query-only MSA; two samples per primary ligand condition",
)]


def sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))


def smiles(path: Path) -> str:
    return json.loads(path.read_text())["PropertyTable"]["Properties"][0]["SMILES"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Continue a live EPSPS BioNeMo run through deterministic workflow stages")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--bionemo-output", type=Path, required=True)
    parser.add_argument("--contacts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    full_sequence = sequence(args.inputs / "P05466.fasta")
    run_summary = json.loads((args.bionemo_output / "run_summary.json").read_text())
    contacts = json.loads(args.contacts.read_text())
    registry = TargetRegistry.default()
    entry = registry.get("AT2G45300", "glyphosate")

    structures = []
    record_scores = {item["record_id"]: item for item in run_summary["records"]}
    for record_id, score in record_scores.items():
        structures.append(StructureModel(
            record_id, entry.agi, "nvidia-bionemo-ir:boltz-2", score["mean_plddt"],
            list(entry.required_context), BIONEMO_PROVENANCE,
            str((args.bionemo_output / f"{record_id}.cif").resolve()), "cif",
            {"mean_plddt": score["mean_plddt"], "ptm": score["ptm"], "iptm": score["iptm"]},
        ))

    def poses(prefix: str, ligand_name: str, key: str) -> list[Pose]:
        result = []
        for record_id, data in contacts["records"].items():
            if record_id.startswith(prefix):
                result.append(Pose(
                    f"{record_id}-{ligand_name}", record_id, ligand_name,
                    record_scores[record_id]["iptm"], data[key], BIONEMO_PROVENANCE,
                    str((args.bionemo_output / f"{record_id}.cif").resolve()),
                    metadata={"contact_cutoff_angstrom": contacts["distance_cutoff_angstrom"]},
                ))
        return result

    s3p_poses = poses("epsps-native", "shikimate-3-phosphate", "s3p_contacts_uniprot")
    backend = PrecomputedComplexBackend(entry.agi, structures, {
        "phosphoenolpyruvate": poses("epsps-native", "phosphoenolpyruvate", "primary_ligand_contacts_uniprot"),
        "shikimate-3-phosphate": s3p_poses,
        "glyphosate": poses("epsps-glyphosate", "glyphosate", "primary_ligand_contacts_uniprot"),
    })
    reasoner = MockScientificBackend()
    workflow = WorkflowOrchestrator(
        registry, backend, backend, backend, backend, reasoner, ArtifactStore(args.output)
    )
    ligand_paths = {name: args.inputs / filename for name, filename in {
        "glyphosate": "glyphosate.json", "phosphoenolpyruvate": "pep.json",
        "shikimate-3-phosphate": "s3p.json",
    }.items()}
    protected = {
        99, 100, 104, 177, 207, 254, 255, 256, 282, 407, 434, 438, 480, 505,
    }
    request = WorkflowRequest(
        target=TargetProtein(entry.agi, entry.protein_name, full_sequence, 1, UNIPROT_PROVENANCE),
        herbicide=Ligand("glyphosate", "herbicide", smiles(ligand_paths["glyphosate"]), provenance=PUBCHEM_PROVENANCE),
        native_ligands=[
            Ligand("phosphoenolpyruvate", "native", smiles(ligand_paths["phosphoenolpyruvate"]), provenance=PUBCHEM_PROVENANCE),
            Ligand("shikimate-3-phosphate", "native", smiles(ligand_paths["shikimate-3-phosphate"]), provenance=PUBCHEM_PROVENANCE),
        ],
        functional_context=[FunctionalPartner(name, "required-context", UNIPROT_PROVENANCE) for name in entry.required_context],
        protected_residues=protected,
    )
    result = workflow.run(request)
    summary = {
        "target": entry.agi, "protein": entry.protein_name, "sequence_source": "UniProt P05466",
        "structure_models": len(result.structures), "pose_hypotheses": len(result.poses),
        "candidate_count": len(result.packets),
        "pareto_front_1": [item.mutation for item in result.pareto_ranking if item.front == 1],
        "candidate_statuses": sorted({packet.status for packet in result.packets}),
        "review_backend": "synthetic fixture; live GPT-Rosalind unavailable",
        "learning_status": "not run: no experimental assay results supplied",
        "limitations": [
            "Query-only MSA", "Boltz-2 poses are hypotheses, not DiffDock or free-energy estimates",
            "No homolog alignment conservation scores", "No Rosetta or molecular-dynamics fold estimates",
            "No live GPT-Rosalind credentials/transport", "No assay data for Milestone 5 recalibration",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "live_run_summary.json").write_text(json.dumps(summary, indent=2))
    (args.output / "protected_residues.json").write_text(json.dumps({
        "residues": sorted(protected), "source": "UniProt P05466 active/binding-site annotations",
    }, indent=2))
    (args.output / "input_manifest.json").write_text(json.dumps({
        "target": {"agi": entry.agi, "uniprot": "P05466", "sequence_length": len(full_sequence),
                   "sequence_sha256": hashlib.sha256(full_sequence.encode()).hexdigest(),
                   "mature_chain": {"start": 77, "end": 520, "structure_numbering_offset": 76}},
        "ligands": {
            "glyphosate": {"pubchem_cid": 3496, "smiles": smiles(ligand_paths["glyphosate"])},
            "phosphoenolpyruvate": {"pubchem_cid": 1005, "smiles": smiles(ligand_paths["phosphoenolpyruvate"])},
            "shikimate-3-phosphate": {"pubchem_cid": 121947, "smiles": smiles(ligand_paths["shikimate-3-phosphate"])},
        },
        "sources": ["https://rest.uniprot.org/uniprotkb/P05466", "https://pubchem.ncbi.nlm.nih.gov"],
        "model": {"name": "boltz-2", "sampling_steps": 100, "samples_per_condition": 2,
                  "msa": "query-only", "contact_cutoff_angstrom": contacts["distance_cutoff_angstrom"]},
    }, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
