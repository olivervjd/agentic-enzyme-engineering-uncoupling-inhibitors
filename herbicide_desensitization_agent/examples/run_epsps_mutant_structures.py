"""Run contact-selected EPSPS substitutions using the installed BioNeMo IR Boltz-2 runtime.

This predicts structures only. No affinity is synthesized from confidence scores.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from herbicide_desensitization_agent.app.validators.mutation_validator import validate_mutation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=2)
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("At least two predictions per ligand condition are required")
    if args.output.exists():
        parser.error("Use a new output directory to avoid mixing predictions")
    from bionemo_ir.data.schemas import InputRequest, MSARecord, Polymer
    from bionemo_ir.pipeline.processor.engine_proc import EngineProcessorConfig, build_processor
    from bionemo_ir.pipeline.stages.configs import FeatureGeneratorStageConfig, WriterStageConfig

    sequence = "".join(line.strip() for line in (args.inputs / "P05466.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    candidates = json.loads(args.candidates.read_text())
    for candidate in candidates:
        if candidate["agi"] != "AT2G45300" or candidate["classification"] != "HERBICIDE_SELECTIVE_MUTABLE":
            raise ValueError("Only EPSPS herbicide-contact single substitutions are accepted")
        validate_mutation(candidate["mutation"], sequence, set())
    mutations = ["WT"] + [item["mutation"] for item in candidates]
    if len(set(mutations)) != len(mutations):
        raise ValueError("Duplicate mutation")
    ligands = {name: json.loads((args.inputs / filename).read_text())["PropertyTable"]["Properties"][0]["SMILES"]
               for name, filename in {"pep": "pep.json", "s3p": "s3p.json", "glyphosate": "glyphosate.json"}.items()}
    args.output.mkdir(parents=True)
    manifest = {
        "target_agi": "AT2G45300", "sequence": sequence,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "modeled_domain": [77, 520], "structure_offset": 76, "chain": "A",
        "model": "boltz-2", "sampling_steps": 100, "msa": "query-only",
        "legend": "Structure predictions for WT and contact-selected EPSPS mutants in PEP+S3P and glyphosate+S3P conditions. "
                  "Replicates use separately seeded runs; these are structural hypotheses, not binding affinities or validated resistance.",
        "records": [],
    }
    for replicate in range(1, args.replicates + 1):
        requests = []
        for mutation in mutations:
            full = sequence
            if mutation != "WT":
                _, position, destination = validate_mutation(mutation, sequence, set())
                if position < 77:
                    raise ValueError("Mutation lies outside the modeled mature chain")
                full = full[:position - 1] + destination + full[position:]
            mature = full[76:]
            for condition, ligand in (("native", "pep"), ("herbicide", "glyphosate")):
                record_id = f"{mutation}-{condition}-{replicate}"
                request = InputRequest(input_id=record_id, polymers=[
                    Polymer(polymer_type="protein", chain_id=["A"], sequence=mature,
                            msas=[MSARecord(content=f">{record_id}\n{mature}\n")]),
                    Polymer(polymer_type="smiles_ligand", chain_id=["B"], sequence=ligands[ligand]),
                    Polymer(polymer_type="smiles_ligand", chain_id=["C"], sequence=ligands["s3p"]),
                ])
                requests.append({"record": request, "__record_id": record_id})
                manifest["records"].append({"record_id": record_id, "mutation": mutation, "condition": condition,
                                             "replicate": replicate, "seed": 42 + replicate,
                                             "sequence_sha256": hashlib.sha256(full.encode()).hexdigest()})
        (args.output / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2))
        processor = build_processor(EngineProcessorConfig(
            model_source="boltz-2", runtime_args={"num_sampling_steps": 100},
            feature_generator_stage=FeatureGeneratorStageConfig(init_context={"random_seed": 42 + replicate}),
            writer_stage=WriterStageConfig(output_path=str(args.output), format="cif"),
            engine_kwargs={"profile_inference": True},
        ))
        rows = list(processor(requests))
        expected = {row["__record_id"] for row in requests}
        if {row.get("__record_id") for row in rows} != expected:
            raise RuntimeError("Inference did not return every requested model")
        for row in rows:
            error = row.get("__inference_error__")
            if isinstance(error, str):
                error = json.loads(error)
            if error and error.get("error_msg"):
                raise RuntimeError(f"Inference failed for {row['__record_id']}: {error['error_msg']}")
            if not (args.output / f"{row['__record_id']}.cif").is_file():
                raise RuntimeError("Missing predicted coordinates")
        print(f"Completed replicate {replicate}: {len(rows)} structures", flush=True)
    print(json.dumps({"output": str(args.output), "structure_count": len(manifest["records"])}))


if __name__ == "__main__":
    main()
