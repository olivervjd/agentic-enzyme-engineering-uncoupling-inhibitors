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
