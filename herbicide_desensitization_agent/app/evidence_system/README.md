# Evidence-system import contract

This package has no scientific-model, network or language-model dependency. It
validates imported tool/experimental measurements, derives contact geometry,
and applies deterministic decisions. Unknown values are `null`; zero represents
an observed zero only. Every API returns JSON-serializable dictionaries.

## Public APIs

- `compare_contexts(left, right, mutant=False, cross_method=False)` checks all
  chemical, biological, mapping, MSA, model, seed and sampling fields. Mutant
  mode requires explicit substitutions matching the canonical sequences.
  Cross-method mode permits different method identities and seeds, but never
  unknown chemistry, sampling differences or changed ligand states.
- `calibrate_thresholds(spec, wt_replicates, controls, registered_at=...,
  candidates_started_at=...)` derives WT envelopes, clips them against the
  preregistered scientific acceptance limits, validates controls and freezes a
  digest. Binding WT ranges must be covered by the declared effect/noise bound.
  Control labels alone are insufficient: resistant/catalytic-loss controls need
  `ligand_changes`, while positive/negative controls need `ligand_affinities`.
  Both are per-evaluation `{value, metric, unit}` measurements with the control's
  scientific provenance. The digest detects edits; revalidation reconstructs
  all derived values from the saved specification and input measurements.
- `assess_binding(record, calibration)` implements the seven evidence states.
  Each independent method needs replicated artifacts, matched contexts,
  physical affinity with a documented interval, pose/contact reproducibility
  and control comparisons. Weakening additionally needs WT-relative evidence
  exceeding the calibrated threshold throughout its interval.
- `affinity_difference(wt, mutant)` subtracts matching endpoint intervals
  conservatively. It never converts between units or physical quantities.
- `evaluate_candidate(candidate, calibration)` returns `{decision, gates,
  limitations, claim}`. The frozen registry controls required ligand/complex
  evaluations. All native intervals, structural metrics, fold stability,
  contact/context/protection evidence and control-derived thresholds must pass.
  The evidence assessment must match the exact variant, ligand, context and
  WT/mutant affinity record; unrelated evidence cannot justify a screen pass.
- `parse_pdb_atoms`, `extract_pose_contacts` and `summarize_contacts` provide
  explicit 4.5 Å heavy-atom contact hypotheses with canonical mapping, atom,
  chain, pose, replicate, method and provenance records. The denominator must
  include every registered pose, including models with no observed contact.
- `classify_residues` uses all native ligands, uncertainty and cross-method
  reproducibility. Missing native contacts are not zero. Second-shell candidates
  require coordinate-based distance evidence. Protected/shared/native-critical
  and highly conserved residues cannot be nominated by this filter.
- `validate_mutation_rationale` requires a falsifiable interaction hypothesis,
  retained native contacts, structural/evolutionary evidence and a scientifically
  sourced stability estimate. `validate_stage` guards stage boundaries.

`tests/test_evidence_system.py` contains complete **synthetic software-only**
examples of the context, calibration, binding-assessment and candidate schemas.
Those fixture values must never be imported as experimental or computational
results for a real protein.

## Integration and trust boundaries

`WorkflowOrchestrator(require_live_gates=True, evidence_calibration=...)` blocks
mutation generation until `validate_live_registry` passes. It requires a
calibrated baseline for the exact target/herbicide/native registry, documented
baseline evidence and (for EPSPS) the specified multicomponent/apo controls.
Existing `require_live_gates=False` workflows remain synthetic demonstrations.

Provenance requires `source_type` (`computational` or `experimental`), `source`
and `artifact`. Artifact references and method-independence families are supplied
by trusted importers. These dictionary validators do not authenticate external
servers, verify arbitrary remote artifact bytes, establish true statistical
independence, or make a self-signed calibration digest a trusted scientific
approval. Archive importers should verify artifact hashes and retain source
files; calibration registration belongs in an access-controlled run registry.
No language-model numeric source can pass a scientific gate.

The reported affinity interval is supplied by the validated predictor/assay.
Subtraction of marginal intervals does not assert joint confidence coverage.
Contact Wilson intervals assume independence of the registered model units.
Uncertain model independence, uncalibrated predictors, pIC50 assay dependence,
homolog-to-target transfer and missing biology remain limitations to disclose.

Unmodified PDB coordinates lack sufficient chemical typing for unambiguous
hydrogen bonds, aromatic contacts or water bridges. Hydrogen bonds require
typed donors/acceptors and explicit D–H–A geometry; salts and hydrophobic
contacts require imported atom typing. Metal contacts and clashes are geometric
hypotheses only. See `CONTACT_DEFINITIONS` for cutoffs and unavailable classes.

A computational screen pass is not experimental herbicide resistance. Reviewed,
provenance-bearing assays covering herbicide interference **and** native enzyme
function are necessary for `EXPERIMENTALLY_VALIDATED`; simply loading an assay
does not approve a candidate.
