# Control-First EPSPS Workflow

## Architecture

`run_integrated_epsps` invokes `WorkflowOrchestrator.run`, not a separate
mutation-screen orchestrator. The calibrated backend can launch Boltz or read
a completed, hash-validated calibration campaign. The original six-target
synthetic entrypoints remain available for regression testing.

| Component | Live implementation |
|---|---|
| Orchestrator | `WorkflowOrchestrator`, with stage manifests and pre-mutation gates |
| Evidence | Bounded Europe PMC abstract retrieval, hashed source records, cited synthesis with explicit homolog applicability |
| Structure quality | `StructureQualityAgent`: MSA, chemical context, WT variability, accuracy and control readiness |
| Pose ensemble | Boltz cofolded substrate complexes; optional public DiffDock-L diagnostic via `PoseEnsembleAgent` |
| Fingerprint | Existing residue-contact comparison, preserving upstream provenance |
| Mutation | Existing conservative single substitutions, reached only after quality and pocket gates |
| Multi-oracle scoring | Existing explicitly heuristic score components, plus protocol-matched candidate measurements; pose confidence is not affinity |
| Review packet | Actual score evidence, mutant measurements and retention-gate results; no hardcoded synthetic labels for real evidence |
| Learning | Integrated assay workflow via `--assays`; pending real assay data, not simulated automatically |
| Judge | Separate model transport from reviewer, plus deterministic evaluation |

Table legend: implementation is distinct from successful scientific validation
or live model access. A stage that fails a gate does not run downstream design.
Unknown evidence is not a pass. No model can authorize assays or bypass a gate.
Quality or evidence failure now permits valid-input pose/contact diagnostics and
workflow-level review/judging to finish. It still blocks mutation generation and
candidate advancement. Workflow reviews are distinct from candidate reviews.

## Calibration Protocol

The calibration runner uses three seeds (101, 103, 107) and upstream Boltz 2.2.1:

- Arabidopsis WT, P05466 residues 77-520: PEP+S3P, glyphosate+S3P, and S3P affinity
  in the PEP+S3P complex.
- E. coli WT, observed protein sequence from experimental PDB 1G6S: PEP+S3P and
  glyphosate+S3P. Only the latter has a matched experimental structural reference.
- E. coli G96A, observed sequence from PDB 1MI4: PEP+S3P, glyphosate+S3P and
  S3P alone. Only S3P alone matches this crystal's ligand context.

This yields 24 fresh structure/affinity predictions, 24 within-context
repeatability comparisons and six experimental-control comparisons. No
experimental structure is supplied as a template. Reference numbering is
explicitly mapped to the observed sequence, and the control sequence difference
is checked. These homolog controls are not substituted into the Arabidopsis
target registry or interpreted as evidence that an Arabidopsis mutation works.

The ColabFold MMseqs2 service receives the two public WT protein sequences.
Returned A3M files are saved, hashed and checked for exact query sequence,
alignment width and distinct homologs. Query-only fallback is prohibited.
Identical MSAs are reused across seeded runs of a sequence. This tests sampling
variability, not evolutionary-input uncertainty. The public service's database
version is not pinned; the downloaded alignment bytes are retained for replay.
The E. coli mutant uses the WT alignment with only the query sequence changed,
keeping homolog rows identical. Historical experimental controls are not
verified as held out from Boltz's training data, so structural agreement is a
sanity check rather than a held-out generalization result.

RDKit audits ligand identity, formula, charge and stereochemistry. This does not
establish the correct protonation state at an unspecified assay pH. The initial
calibration deliberately retains the recorded chemical states to avoid mixing
an unvalidated chemical change into the MSA benchmark. Chemistry readiness
remains false until an appropriate protocol is supplied and independently checked.

Original structural thresholds are reported unchanged. Three-seed ranges are
descriptive, not confidence intervals. Boltz pIC50 remains separate from Kd,
Km, Ki, catalytic activity and the final Kd-based function-retention gate.
Known-control direction checks are diagnostic, not quantitative calibration
between different experimental endpoints.

Sources: [Boltz inputs and MSA](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md),
[1G6S WT reference](https://www.rcsb.org/structure/1G6S),
[1MI4 control and associated study](https://www.rcsb.org/structure/1MI4).

## Run

Use the existing isolated Linux Boltz environment and model cache. Structural
analysis also requires the `calibration` optional dependency set. Do not install
Boltz into the separate BioNeMo environment.

```bash
python -m herbicide_desensitization_agent.examples.run_integrated_epsps \
  --previous /path/to/completed-docking-affinity-campaign \
  --calibration /path/to/new-calibration \
  --output /path/to/new-agent-workflow \
  --boltz-cache /path/to/boltz-cache \
  --boltz /path/to/boltz-environment/bin/boltz \
  --diffdock-cache /path/to/diffdock-public-cache
```

A completed calibration directory can be reused only when its output artifacts
match the recorded hashes. The workflow output directory must be fresh.
`run_epsps_calibration` exposes the same calibration backend independently;
its `--analyze-only` option rebuilds reports after computation without changing
input requests or rerunning predictions.
For correcting the control MSA in an earlier calibration, `run_epsps_calibration
--reuse-calibration /path/to/earlier-calibration` accepts a fresh output directory,
validates the source artifacts, requests, MSA bytes, seeds and inference settings,
copies the 15 unchanged WT predictions, and reruns the nine mutant-control
predictions with matched homolog rows. Reused and new counts are reported separately.

Exit code 2 means `NEEDS_REASSESSMENT`, not an infrastructure crash or successful
candidate screen. `execution.json` and `workflow_stages.json` record the reason.
Other exceptions are failures. A successful computer calculation is not a
scientific pass.

## Separate Models

Supply `OPENAI_API_KEY` securely in the process environment. Set `EVIDENCE_MODEL`,
`REVIEW_MODEL` and `JUDGE_MODEL` to exact API model IDs that the account can access,
or use the corresponding CLI arguments. No model ID is guessed from a Codex UI
label and no model is silently substituted. The same key may authenticate these
transports, but each role has its own configured model. Responses use strict
structured output and `store: false`; model requests still transmit the supplied
evidence to OpenAI. No credential is written to artifacts.

The evidence model now defaults to the documented `gpt-rosalind-research` API ID.
Use `--use-codex-api-key` to explicitly reuse an actual Codex API-key login, or
provide `OPENAI_API_KEY`. A missing evidence key is an error, not an implicit
curated fallback. Use `--curated-evidence-only` only when deliberately choosing
no model synthesis. Reviewer and judge default to `gpt-5.6-luna`, with separate
requests and configurable model IDs. An access preflight records unavailable
models without blocking independent diagnostics or silently substituting models.
Cited curated evidence and deterministic packets remain distinguishable from
live LLM synthesis/review. In live mode, an absent judge never falls back to the
reviewer. Configured calls that fail or are refused raise an error, rather than
generating mock responses. Model availability must be verified with actual
account access; transport fixture tests do not establish that access.

The September 20 access-only preflight returned HTTP 404 for
`gpt-rosalind-research` using the existing Codex API key. It was not listed in that
key's model catalog. The default is configured, but Rosalind inference is **not
enabled by the account**. The earlier check used `gpt-rosalind`, not the documented
API ID. Current [OpenAI documentation](https://developers.openai.com/api/docs/pricing)
limits Rosalind to approved internal research through trusted access.

The evidence adapter now retrieves up to eight relevance-ranked Europe PMC
records with abstracts before synthesis, saving the query, timestamp, response
bytes, source hashes and applicability limitations. Retrieval does not require
Rosalind access. The requested Rosalind synthesis remains unavailable until
approved; retrieved text is untrusted data, not instructions. This is not an
exhaustive literature review, full-text extraction or automatic binding-label
curation. Use `--no-literature-retrieval` to explicitly disable online retrieval.

For a bounded check using existing predictions, add `--diagnostics-only
--use-codex-api-key`. This never launches Boltz or generates mutations. Optional
`--diffdock-cache` still requests independent diagnostic docking; omit it to
avoid new docking jobs. `--cached-only` also prevents Boltz launches but allows
design after all scientific gates pass. Diagnostic reviews use two real model
requests when available and cannot approve candidates.

`--assays /path/to/assays.json` accepts a JSON list matching `AssayResult`, with
measurement units and provenance. Learning validates the target/mutation against
reviewed candidates before recalibration. Missing assays produce `PENDING_ASSAYS`,
not a completed training claim; explicitly synthetic provenance is rejected in
live mode. `workflow_review_packet.json`, `workflow_review.json` and
`learning_readiness.json` document the result even when design is blocked.

Access-only check, with no workflow or inference run:

```bash
python -m herbicide_desensitization_agent.examples.check_model_access \
  --use-codex-api-key --judge-model gpt-5.6-luna \
  --output /path/to/model-access-preflight.json
```

## Remaining Scientific Gates

The current calibration cannot authorize new candidates solely from these
outputs. It lacks a validated protonation protocol, a same-target experimental
accuracy benchmark, quantitatively appropriate functional control calibration,
and independent S3P-aware pocket confirmation. WT repeatability is calculated
from the new outputs and may also fail.

Public DiffDock-L cannot represent the fixed S3P cosubstrate in this workflow.
Its adapter deliberately strips non-protein receptor atoms and labels every
result diagnostic-only. Adding S3P to a file that the model ignores is not a
fix. The pocket gate must remain closed until a validated context-capable
independent method supplies evidence. Optional diagnostics can now run despite
an unresolved structure-quality gate; they never override a failing gate.

`ArtifactCandidateOracles` requires mutation-specific, protocol-matched artifacts
and measurements when the earlier gates pass. It does not reuse the previous
query-only mutation campaign as if it were calibrated. Automatic regeneration
of a complete new mutant campaign and a calibrated affinity/stability oracle
are still required before claiming an end-to-end autonomous design loop.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests cover source integrity, citation grounding, valid MSAs, sequence binding,
quality/pocket blocking, model-role separation, actual review evidence, missing
mutant oracles, structured API payloads, refusals and sanitized API errors.
