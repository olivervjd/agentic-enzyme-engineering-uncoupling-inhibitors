# Calibrated evidence workbench

The evidence workbench connects raw molecular models, per-ligand affinity outputs,
residue contacts, control calibration, candidate decision gates and the next
experiment. Every live mutation-nomination run now requires a prospectively
frozen calibration in addition to the original structure-quality checks.
Missing calibration blocks design; it does not trigger threshold relaxation.

The recorded EPSPS example contains 24 matched-MSA WT/homolog-control predictions
and a separately labelled archive of 42 older WT/mutant predictions. It does not
establish an experimentally validated mutant. Original ligand quantities remain
separate, and absent prediction intervals are never replaced by seed ranges.

## Build from recorded evidence

Use Python 3.11 or later and the optional scientific dependencies:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[calibration,reports]'
.venv/bin/python -m herbicide_desensitization_agent.examples.export_evidence_dashboard \
  --calibration /path/to/epsps-calibration-matchedmsa-20260920 \
  --archive /path/to/full-workflow-20260919-docking-affinity \
  --workflow /path/to/epsps-integrated-matchedmsa-20260920 \
  --output /path/to/evidence-workbench
```

The exporter verifies recorded calibration hashes and all archived coordinate
hashes, recalculates atom-level proximity contacts, conservatively retains
missing evidence, and writes five CSV/JSON tables with legends. It computes
replicate pocket side-chain RMSD, ligand displacement under matched atom names,
homolog experimental contact recovery, and descriptive MSA query-identity
fractions. The geometric definitions are recorded with results. Unsupported
chemical interaction types, ligand atom mappings, pocket-volume and stability
estimates stay missing.

Copy the exported `data/` directory into `demo/public/data/`, then:

```sh
cd demo
npm ci
npm run build
cd ..
.venv/bin/python -m herbicide_desensitization_agent.examples.package_evidence_dashboard \
  --output /path/to/evidence-workbench --dist demo/dist
```

`report.html` is a self-contained interactive report with embedded data, molecular
geometry and runtime libraries. It can be opened offline. `dashboard/` can be
served by a local HTTP server. Raw evidence remains separately downloadable in
the bundle, indexed by `workflow-manifest.json`. Figures are 300 dpi PNG and
editable vector SVG. No runtime authentication material enters either report.

## Prospective workflow

See [evidence contracts](evidence_contracts.md) for deterministic calibration,
chemical-context, binding, residue and candidate schemas. Pass a reviewed frozen
record with `run_integrated_epsps --evidence-calibration calibration.json`.
The live orchestrator verifies registry scope and complete EPSPS multicomponent
contexts before nomination. Calibration is derived from WT/control records;
an existing candidate result must never be used to retrospectively loosen it.
The five other fixed registry pairings require their own validated baseline;
the EPSPS data do not fill those gaps.

The new evidence states are separate from legacy function-retention results.
A computational pass means `MEETS_COMPUTATIONAL_SCREEN`, never herbicide
resistance. `EXPERIMENTALLY_VALIDATED` requires actual provenance-bearing
measurements; an LLM opinion cannot establish it.

## Explanatory model review

Evidence, reviewer and judge defaults are `gpt-5.6-luna`. The existing Codex API
key is reused only by explicit opt-in, in memory. The exporter is deterministic
and does not require an API key. Optional model commentary is generated separately:

```sh
.venv/bin/python -m herbicide_desensitization_agent.examples.review_evidence_bundle \
  --input /path/to/evidence-workbench/data/evidence_system.json --use-codex-api-key
```

This writes a separate redacted review record. Neither reviewer nor judge can
modify scientific fields or approve candidates. A changed source file marks the
review stale. Metadata access, inference success, and scientific validation are
different statuses. See the [official Luna model documentation](https://developers.openai.com/api/docs/models/gpt-5.6-luna)
and [Codex authentication documentation](https://learn.chatgpt.com/docs/auth).

## Current scientific limitations

Current Arabidopsis results lack a validated assay pH/protonation protocol,
same-target control panel, matched-context independent binding evidence,
calibrated affinity intervals, fold stability and experimental native-function
measurements. Native and herbicide contact sets overlap. These gaps keep the
nomination gate closed. Apo controls, control-panel calibration, independent
redocking/cross-docking, shortlist free-energy/MD and eventual biochemical
measurements are represented as explicit pending work, not fabricated results.
No legal gene-editing eligibility is claimed without a specified jurisdiction,
organism, nucleotide edit and applicable framework.

## Verification

```sh
.venv/bin/python -m unittest discover -s tests -q
cd demo
npm run build
npm run test:browser
```

Synthetic positive and negative tests exercise the contracts; they are software
fixtures and never appear in the scientific dashboard.
