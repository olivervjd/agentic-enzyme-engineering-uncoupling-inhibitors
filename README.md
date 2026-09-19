# Herbicide Desensitization Agent

Milestone 1 prototype for a review-gated computational decision engine covering
six fixed *Arabidopsis thaliana* AGI–herbicide pairings across five herbicide
classes.

The current implementation contains a fixed provenance-aware target registry,
typed schemas, deterministic validators, an end-to-end orchestrator, and six
placeholder examples. Milestone 2 adds a BioNeMo IR `build_processor` adapter,
a DiffDock NIM HTTP adapter, ensemble manifests, geometric contact extraction,
interaction fingerprints, and target-specific structure-context gates.
Each persisted BioNeMo structure also produces a standalone Mol* HTML viewer,
a copied coordinate file, a score manifest, and `visualizations.json`.

It does **not** predict experimentally actionable mutations. Mock outputs are
synthetic plumbing fixtures and must not be treated as biological evidence.

## Supported targets

| Herbicide | Protein | AGI |
|---|---|---|
| Glyphosate | EPSPS | AT2G45300 |
| Glufosinate | GLN1;2 | AT1G66200 |
| Glufosinate | GLN2 | AT5G35630 |
| Atrazine | PsbA | ATCG00020 |
| Chlorsulfuron | ALS/AHAS (CSR1) | AT3G48560 |
| 2,4-D | TIR1 | AT3G62980 |

## Run

```bash
python -m herbicide_desensitization_agent.examples.run_all
python -m unittest discover -s tests -v
python -m herbicide_desensitization_agent.examples.render_milestone2_validation \
  output/live-smoke/0.cif output/live-smoke/0_scores.json --output artifacts
```

Real NVIDIA BioNeMo/NIM and GPT-Rosalind clients should be added behind the
interfaces in `app/backends/`; the orchestrator does not depend on vendor SDKs.

## Milestone 2 integration

`BioNeMoIRStructureBackend.from_runtime(...)` uses NVIDIA's supported BioNeMo
Inference Runtime processor surface. It requires the optional `bionemo-ir`
package, compatible NVIDIA hardware, model weights, and an explicitly prepared
GPU environment. OpenFold2, OpenFold3, and Boltz model keys are configuration,
not hard-coded dependencies.

`DiffDockNIMBackend` calls the deployed NIM endpoint
`/molecular-docking/diffdock/generate` through an injected JSON transport. Check
the deployment's `/docs` schema and provide a response adapter if its response
envelope differs. Credentials are supplied at runtime and are never persisted by
this repository.

Monomer-only BioNeMo requests intentionally fail the PsbA and TIR1 biological-
context gates. Callers must supply complex-aware request factories and declare
only context actually represented in the generated model.

Mol* viewer HTML embeds the predicted coordinates so it can be opened directly.
The Mol* JavaScript and CSS are loaded from jsDelivr, so an internet connection
is required when viewing the file. These viewers display model output and
confidence metadata; they do not establish biological validity.

## Milestone 3 mutation and scoring

Milestone 3 generates conservative single substitutions at fingerprint-supported
herbicide-selective contacts by default. Explicitly enabling second-shell selection
requires spatial evidence in pose metadata; sequence adjacency is not spatial evidence. Shared, native-critical,
and explicitly protected positions are excluded. Nine score components retain
their own evidence and uncertainty, and candidates are assigned non-dominated
Pareto fronts without collapsing them to one opaque score.

The built-in conservation and fold scores accept homolog substitution frequency
and predicted delta-delta-G metadata when external tools provide them. Otherwise
they emit explicit neutral/physicochemical baselines and uncertainty warnings.
Generate all fixed-pairing validation artifacts with:

```bash
python -m herbicide_desensitization_agent.examples.run_milestone3_validation \
  --output artifacts/milestone-3
```

## Milestone 4 review and evaluation

Milestone 4 expands every candidate into a structured review packet and evaluates
it with a GPT-Rosalind-compatible twelve-category domain judge plus an independent
deterministic judge. The model judge is transport-injected; credentials and model
responses are not persisted by the adapter. Deterministic validation runs before
judging and prevents automatic assay approval.

The fixed-target benchmark runner rejects out-of-registry pairings and reports
the eight required retrospective metrics. Bundled labels are explicitly synthetic;
real recall claims require a separately curated, provenance-preserving hidden set.
Ten adversarial conditions are mapped to required safe handling.

```bash
python -m herbicide_desensitization_agent.examples.run_milestone4_validation \
  --output artifacts/milestone-4
```

## Milestone 5 learning loop

Milestone 5 ingests provenance-bearing assay measurements, enforces mechanism-
specific native-function readouts, estimates component-wise prediction residuals,
and selects unassayed single substitutions using a transparent exploitation/
uncertainty acquisition rule. It never changes a candidate to an approved state.

The recalibration is a bounded residual correction, not a causal biological model.
Small datasets carry an explicit instability warning. Required outputs are
`assay_results.json`, `model_recalibration_report.json`, and
`next_round_candidates.json`.

```bash
python -m herbicide_desensitization_agent.examples.run_milestone5_validation \
  --output artifacts/milestone-5
```

## Continuing a live EPSPS run

The precomputed-complex adapter can continue provenance-bearing BioNeMo outputs
through contact fingerprinting, mutation generation, scoring, Pareto ranking,
visualization, and governed review without rerunning inference:

```bash
python -m herbicide_desensitization_agent.examples.run_epsps_precomputed_live \
  --inputs work/epsps-real-inputs \
  --bionemo-output outputs/epsps-live-run/bionemo/epsps-live-output \
  --contacts outputs/epsps-live-run/contact_report.json \
  --output outputs/epsps-live-run/workflow
```

This command does not substitute mocks for unavailable scientific services. A
run using the mock reasoning fixture is labeled as such, and the learning stage
remains pending until real assay data are provided.

## Function-retention negative-design gate

Every collection receives `function_retention_report.json`, CSV, and Markdown
tables with a WT reference row, units, prediction intervals, methods, and a table
legend. A companion legend and threshold manifest are saved beside the CSV.
The gate requires mutant-versus-reference TM-align similarity, CA lDDT and
coverage, global CA and active-site backbone RMSD, direct Kd estimates for
herbicide and every registered native ligand, prediction intervals, matching
affinity protocols, and fold ΔΔG. WT-relative ratios are calculated from the
actual estimates; inconsistent supplied ratios and non-finite values are rejected.
Heuristic contact scores and Boltz pIC50 are never relabeled or converted to Kd.

The default screen requires both directional TM-scores ≥0.95, CA lDDT ≥0.90,
coverage ≥0.95, global CA and active-site backbone RMSD ≤1.0 Å, and fold ΔΔG
≤2.0 kcal/mol. The entire herbicide Kd ratio interval must be ≥10; every native
ligand ratio interval must be inside [1/3, 3], excluding both excessive weakening
and excessive strengthening. Ratio intervals conservatively combine supplied
bounds, not just overlapping point estimates. Thresholds are configurable triage
criteria, not a claim of perfectly identical structure or absent binding.
Known failures produce `FAILS_COMPUTATIONAL_SCREEN` even when other evidence is
missing. Otherwise missing evidence produces `INSUFFICIENT_EVIDENCE`, never a
pass. These decisions now gate the review packets as well as the tables.

### Coordinate comparisons and evidence

Install the optional structural dependencies with `python -m pip install '.[structure]'`.
`TMAlignStructuralMatcher` uses the existing TM-align implementation in `tmtools`
and Biopython's Kabsch fit. Chain IDs, numbering offsets and modeled sequence
domains are explicit. The expected single mutation is verified against the
coordinate sequence; WT coordinates cannot masquerade as mutant predictions.
Active-site N/CA/C/O RMSD uses the global CA transform without a separate site fit.
Missing site atoms yield missing evidence, not zero RMSD. All provided comparisons
are aggregated conservatively, taking minimum similarity and maximum RMSD.

Pass `--retention-evidence evidence.json` to `run_epsps_precomputed_live` to use
real evidence. The JSON is keyed by `WT` and mutation identity. Each record uses:

- `target_agi`, `provenance` (source, method, evidence_type), and `structural_method`.
- `ligand_kd_molar`: ligand name to positive Kd in M.
- `ligand_kd_intervals_molar`: ligand name to `[lower, upper]` in M.
- `affinity_protocol`: ligand name to an identical protocol identifier for WT and
  mutant, including model/version, chemical state, conditions, and uncertainty method.
- `interval_description`: how the supplied bounds were obtained and what they mean.
- `fold_ddg_kcal_mol` and either `structural_metrics` or `structure_comparisons`.

Each `structure_comparisons` item contains `mutant`, `reference`, `mutant_chain`,
`reference_chain`, `mutant_offset`, `reference_offset`, `active_site_residues`,
`sequence_start`, and `sequence_end`. Paths resolve relative to the evidence file.
For the existing EPSPS mature-chain models, use chain A, offsets 76, and sequence
domain 77-520. Residue selections use full sequence numbering. WT comparisons
should use the same reference on both sides as an explicitly labeled identity
control. Compare WT replicates separately to assess model variability.

### EPSPS mutant structure run

The installed H200 BioNeMo IR environment can run the contact-selected candidates:

```bash
python -m herbicide_desensitization_agent.examples.run_epsps_mutant_structures \
  --inputs work/epsps-real-inputs \
  --candidates work/contact-only-workflow/AT2G45300-glyphosate/candidate_mutations.json \
  --output work/epsps-mutant-structures --replicates 2
```

This creates WT and mutant structures in both PEP+S3P and glyphosate+S3P
conditions, with separately seeded replicate runs and a sequence manifest. It
does not calculate Kd or folding ΔΔG. Query-only MSA predictions are preliminary.
The repository currently has real starting complexes for EPSPS only; other
target examples remain synthetic. Native-binding checks cover the registry's
listed ligands, not every aspect of enzyme turnover, cofactors, assembly, or
TIR1 signaling. Additional validated evidence is required to establish those.

Methods: [TM-align implementation](https://github.com/jvkersch/tmtools),
[Foldseek alignment options](https://github.com/steineggerlab/foldseek), and
[Boltz affinity interpretation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).
