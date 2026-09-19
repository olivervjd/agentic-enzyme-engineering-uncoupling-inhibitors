"""Validate completed EPSPS predictions and report every comparison, including WT controls."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
from itertools import combinations
import json
from pathlib import Path

from herbicide_desensitization_agent.app.agents.function_retention import (
    FunctionRetentionAgent, RETENTION_LEGEND, retention_csv, retention_markdown,
)
from herbicide_desensitization_agent.app.backends.precomputed import PrecomputedFunctionRetentionBackend
from herbicide_desensitization_agent.app.backends.structural import TMAlignStructuralMatcher
from herbicide_desensitization_agent.app.schemas.models import (
    Ligand, MutationCandidate, Provenance, ScorePacket, TargetProtein,
)
from herbicide_desensitization_agent.app.storage.artifact_store import ArtifactStore
from herbicide_desensitization_agent.app.validators.mutation_validator import validate_mutation


COMPARISON_LEGEND = (
    "Table/Figure legend: Every row compares two independently predicted EPSPS mature chains "
    "(full sequence residues 77-520; chain A; coordinate numbering offset +76). For each mutation, "
    "all mutant-versus-WT replicate pairs are evaluated within each ligand condition. Native means "
    "PEP+S3P; herbicide means glyphosate+S3P. WT-control rows compare the two WT replicates within "
    "each condition, not WT to itself. The final retention table instead uses a WT identity reference. "
    "TM-align scores are length-normalized structural similarities, not prediction confidence. "
    "CA lDDT measures agreement of mapped CA distances within a 15 A reference neighborhood, using "
    "0.5/1/2/4 A tolerances and averaging per residue. RMSDs use one global Kabsch CA fit; active-site "
    "N/CA/C/O RMSD covers the union of WT ligand-contact positions from the original 5 A contact "
    "report. Lower RMSD and higher TM/lDDT indicate closer models. Each point is one comparison, "
    "not an independent biological observation or confidence interval. WT variability provides "
    "context, not grounds to relax thresholds after seeing the results. Query-only MSA Boltz-2 "
    "predictions do not establish native function or herbicide resistance. No affinity is inferred from these structural comparisons."
)


def read_campaign(directory):
    manifest = json.loads((directory / "campaign_manifest.json").read_text())
    if manifest["target_agi"] != "AT2G45300" or manifest["modeled_domain"] != [77, 520]:
        raise ValueError("This evaluator requires the EPSPS mature-chain campaign")
    sequence = manifest["sequence"]
    if hashlib.sha256(sequence.encode()).hexdigest() != manifest["sequence_sha256"]:
        raise ValueError("Campaign sequence hash mismatch")
    records = manifest["records"]
    ids = [item["record_id"] for item in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate record ID")
    mutations = {item["mutation"] for item in records}
    if "WT" not in mutations:
        raise ValueError("WT reference is required")
    groups = {}
    for item in records:
        identity = item["record_id"]
        if Path(identity).name != identity:
            raise ValueError("Record IDs must be simple filenames")
        if not (directory / f"{identity}.cif").is_file():
            raise ValueError(f"Missing coordinates for {identity}")
        mutation = item["mutation"]
        if item["condition"] not in {"native", "herbicide"}:
            raise ValueError("Unknown ligand condition")
        full = sequence
        if mutation != "WT":
            _, position, destination = validate_mutation(mutation, sequence, set())
            full = full[:position - 1] + destination + full[position:]
        if hashlib.sha256(full.encode()).hexdigest() != item["sequence_sha256"]:
            raise ValueError(f"Mutant sequence hash mismatch for {identity}")
        key = (mutation, item["condition"])
        groups.setdefault(key, []).append(item)
    for mutation in mutations:
        for condition in ("native", "herbicide"):
            rows = groups.get((mutation, condition), [])
            if len(rows) < 2 or len({item["replicate"] for item in rows}) != len(rows):
                raise ValueError(f"Missing/duplicate replicate for {mutation}/{condition}")
            if {item["replicate"] for item in rows} != {item["replicate"] for item in groups[("WT", condition)]}:
                raise ValueError("WT and mutant replicate sets must match")
    return manifest, groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--workflow", type=Path, required=True, help="Contact-only workflow directory containing candidates and scores")
    parser.add_argument("--contacts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--affinity-evidence", type=Path, help="Optional WT/mutant direct Kd evidence; never inferred from structure scores")
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    manifest, groups = read_campaign(campaign)
    raw_candidates = json.loads((args.workflow / "candidate_mutations.json").read_text())
    raw_scores = json.loads((args.workflow / "score_packets.json").read_text())
    candidates = [MutationCandidate(**{**item, "provenance": [Provenance(**p) for p in item["provenance"]]})
                  for item in raw_candidates]
    scores = [ScorePacket(**{**item, "provenance": [Provenance(**p) for p in item["provenance"]]}) for item in raw_scores]
    mutations = ["WT"] + [item.mutation for item in candidates]
    if set(mutations) != {key[0] for key in groups}:
        raise ValueError("Campaign mutations differ from the contact-selected workflow")
    if any(item.classification != "HERBICIDE_SELECTIVE_MUTABLE" or item.agi != "AT2G45300" for item in candidates):
        raise ValueError("Campaign requires contact-selected EPSPS candidates")
    contact_report = json.loads(args.contacts.read_text())
    site = sorted({position for record in contact_report["records"].values()
                   for key, positions in record.items() if key.endswith("_contacts_uniprot") for position in positions})
    if not site:
        raise ValueError("An explicit WT ligand-contact active site is required")
    options = dict(mutant_chain="A", reference_chain="A", mutant_offset=76, reference_offset=76,
                   active_site_residues=site, target_sequence=manifest["sequence"], sequence_start=77, sequence_end=520)
    matcher = TMAlignStructuralMatcher()
    evidence = json.loads(args.affinity_evidence.read_text()) if args.affinity_evidence else {}
    comparisons = []
    for mutation in mutations:
        metrics = []
        for condition in ("native", "herbicide"):
            for query in groups[(mutation, condition)]:
                references = groups[("WT", condition)] if mutation != "WT" else [query]
                for reference in references:
                    values = matcher.compare(campaign / f"{query['record_id']}.cif",
                                             campaign / f"{reference['record_id']}.cif", mutation=mutation, **options)
                    metrics.append(values)
                    comparisons.append({"mutation": mutation, "condition": condition,
                                        "query": query["record_id"], "reference": reference["record_id"], **values})
        record = evidence.setdefault(mutation, {})
        if record.get("target_agi", "AT2G45300") != "AT2G45300":
            raise ValueError("Affinity target mismatch")
        record.pop("structure_comparisons", None)
        record.update(target_agi="AT2G45300", structural_method=matcher.method + "; worst over all within-condition replicate pairs",
                      structural_metrics={key: None if any(row[key] is None for row in metrics) else
                                          (max(row[key] for row in metrics) if "rmsd" in key else min(row[key] for row in metrics))
                                          for key in metrics[0]})
        record.setdefault("provenance", []).append({
            "source": str(campaign / "campaign_manifest.json"), "method": matcher.method,
            "evidence_type": "computed-structure-comparison",
            "notes": "No affinity or ddG inferred from coordinates. Modeled domain: 77-520. Active site: " + str(site),
        })
    for condition in ("native", "herbicide"):
        for reference, query in combinations(groups[("WT", condition)], 2):
            values = matcher.compare(campaign / f"{query['record_id']}.cif",
                                     campaign / f"{reference['record_id']}.cif", mutation="WT", **options)
            comparisons.append({"mutation": "WT-control", "condition": condition,
                                "query": query["record_id"], "reference": reference["record_id"], **values})
    provenance = [Provenance(str(campaign / "campaign_manifest.json"), "Boltz-2 query-only MSA", "predicted-structure")]
    agent = FunctionRetentionAgent(PrecomputedFunctionRetentionBackend(evidence))
    records = agent.assess(TargetProtein("AT2G45300", "EPSPS", manifest["sequence"], provenance=provenance),
                           candidates, scores, Ligand("glyphosate", "herbicide", "", provenance=provenance),
                           [Ligand(name, "native", "", provenance=provenance)
                            for name in ("phosphoenolpyruvate", "shikimate-3-phosphate")], [])
    store = ArtifactStore(args.output)
    run = "AT2G45300-glyphosate"
    store.write_manifest(run, "function_retention_report.json", records)
    store.write_manifest(run, "retention_evidence.json", evidence)
    store.write_manifest(run, "structure_comparisons.json", comparisons)
    store.write_manifest(run, "function_retention_thresholds.json", agent.thresholds)
    store.write_manifest(run, "campaign_manifest.json", manifest)
    native = ["phosphoenolpyruvate", "shikimate-3-phosphate"]
    store.write_text(run, "function_retention_table.csv", retention_csv(records, "glyphosate", native))
    store.write_text(run, "function_retention_table.md", retention_markdown(records, "glyphosate", native))
    store.write_text(run, "figure_and_table_legends.md", RETENTION_LEGEND + "\n\n" + COMPARISON_LEGEND + "\n")
    csv_buffer = io.StringIO()
    writer = csv.DictWriter(csv_buffer, fieldnames=list(comparisons[0]))
    writer.writeheader()
    writer.writerows(comparisons)
    store.write_text(run, "structure_comparisons.csv", csv_buffer.getvalue())
    print(json.dumps({"structures": len(manifest["records"]), "comparisons": len(comparisons),
                      "decisions": {row.mutation: row.decision for row in records},
                      "output": str(store.run_directory(run))}, indent=2))


if __name__ == "__main__":
    main()
