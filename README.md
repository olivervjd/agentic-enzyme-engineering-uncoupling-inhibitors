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

Milestone 3 generates conservative single substitutions only at fingerprint-
supported herbicide-selective or second-shell positions. Shared, native-critical,
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

Every mutation now receives `function_retention_report.json`, CSV, and Markdown
tables. The gate requires mutant-versus-reference Foldseek/TM-align similarity,
alignment LDDT and coverage, active-site RMSD, direct Kd estimates for herbicide
and every native ligand, Kd fold changes relative to wild type, and fold ΔΔG.
Heuristic contact scores and Boltz pIC50 are never relabeled or converted to Kd.

The default screen requires both directional TM-scores ≥0.80, alignment LDDT
≥0.70, coverage ≥0.80, active-site RMSD ≤1.5 Å, at least tenfold weaker herbicide
binding, no more than threefold weaker native-ligand binding, and fold ΔΔG ≤2.0
kcal/mol. These are configurable computational triage thresholds, not biological
proof or authorization for an experiment. Missing direct evidence produces
`INSUFFICIENT_EVIDENCE`, never a pass.
