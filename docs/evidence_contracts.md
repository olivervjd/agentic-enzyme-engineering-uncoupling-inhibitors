# Prospective evidence contracts

Import public functions from
`herbicide_desensitization_agent.app.evidence_system`. Inputs and outputs are
plain JSON-compatible dictionaries; no LLM produces numeric scientific values.
`None`/JSON `null` means unknown. An absent native-ligand observation must never
be imported as a zero contact frequency or preserved binding.

All numbers in this document's examples are **synthetic software examples**,
not measurements, predictions, recommended cutoffs, or a calibrated EPSPS
protocol. Full executable synthetic fixtures are in
`tests/test_evidence_system.py`. They are confined to tests.

## Minimal safe usage with incomplete real data

```python
from herbicide_desensitization_agent.app.evidence_system import (
    assess_binding, evaluate_candidate, validate_stage,
)

binding = assess_binding({
    "variant": "WT", "ligand": "glyphosate", "evaluated": True,
    "methods": [],
}, calibration=None)
assert binding["evidence_state"] == "Insufficient evidence"

decision = evaluate_candidate({"variant": "candidate_1"}, calibration=None)
assert decision["decision"] == "INSUFFICIENT_EVIDENCE"
# The UI must display each decision['gates'] reason, never fill missing values.
```

## Shared objects

Scientific provenance has these fields:

```json
{"source_type": "computational", "source": "validated-tool-and-version", "artifact": "path-or-persistent-artifact-id"}
```

`source_type` accepts `computational` or `experimental`. A language-model source
cannot pass scientific validation. Importers must verify artifact identity and
integrity; the core does not download or authenticate arbitrary source files.
Preserve their SHA256 hashes in the workflow's evidence manifest.

An affinity object requires every field below:

```json
{
  "metric": "pIC50", "unit": "dimensionless", "estimate": 5.0,
  "interval": [4.9, 5.1], "interval_method": "synthetic example only",
  "interval_level": 0.95,
  "provenance": {"source_type": "computational", "source": "synthetic fixture", "artifact": "fixture://test-only"}
}
```

Supported physical quantity/unit pairs are `Kd`, `Ki` or `IC50` with `M`, `mM`,
`uM`, `nM` or `pM`; `pIC50`, `pKd` or `pKi` with `dimensionless`; and `delta_G`
with `kcal/mol` or `kJ/mol`. Matching units are required, even where conversion
would be possible. All concentration bounds must be positive. The interval
must contain the finite estimate. Docking confidence is never an affinity.
`affinity_difference(wt, mutant)` reports mutant minus WT; its interval is
conservative endpoint subtraction, not an assertion of joint statistical
coverage. No conversion among concentration, pIC50, ΔG or confidence occurs.

## Complete context object

`compare_contexts(left, right, mutant=False, cross_method=False)` requires:

| Field | JSON type / requirement |
| --- | --- |
| `protein_sequence` | Canonical sequence string |
| `isoform` | Resolved isoform identifier |
| `modeled_residue_range` | `[first, last]`, ordered positive canonical integers |
| `residue_mapping` | Nonempty mapping, e.g. `{"A:12A": 15}` including insertion codes |
| `oligomeric_state` | Resolved assembly description |
| `cellular_context` | Subcellular/membrane description |
| `substrates`, `cofactors`, `metal_ions`, `catalytic_waters`, `partners` | Explicit lists; `[]` means confirmed absence |
| `ligands` | List of chemical-state objects below |
| `pH` | Finite number between 0 and 14 |
| `ionic_assumptions` | Resolved conditions string |
| `structure_model`, `structure_model_version` | Explicit method and version |
| `affinity_model`, `affinity_model_version` | Explicit method and version |
| `msa_source`, `msa_depth` | Documented source and positive integer depth |
| `seed` | Nonnegative integer |
| `sampling` | Nonempty object of sampling settings |
| `substitutions` | Mutant only, e.g. `["C2A"]`; must exactly match sequence changes |

Each ligand requires `id`, `canonical_identifier`, `stereochemistry`,
`protonation`, `tautomer` (resolved strings) and integer `formal_charge`.
`unknown`, `unspecified`, `unavailable`, `pending`, missing or empty values are
incomplete evidence. An apo complex uses `ligands: []` and `complex_type: apo`.

WT/mutant comparison requires matching versions, seed, sampling and chemistry.
Cross-method comparison permits different model names/versions and seeds, but
requires their recorded identities and matching chemical/biological/sampling
context. Non-equivalent contexts cannot directly validate agreement or conflict.

## Calibration input and temporal scope

Call `calibrate_thresholds(spec, wt_replicates, controls,
registered_at="...Z", candidates_started_at="...Z")` **before inspecting
candidate results**. Both timestamps must be timezone-aware; registration must
precede candidate evaluation. The protocol belongs in an access-controlled
prospective registry. A digest alone is not proof of preregistration.

`spec` requires:

| Field | Meaning |
| --- | --- |
| `metrics` | `{metric_name: {direction: "min" or "max", limit: number, margin: nonnegative_number, unit: string}}` |
| `minimum_wt_replicates` | Optional; at least 3, default 3 |
| `required_evaluations` | Frozen `{evaluation_id: role}` registry with herbicide and native substrate roles |
| `assembly_required` | Boolean; true adds an `assembly_change` gate |
| `binding` | The object detailed below |
| `contact` | Optional frozen differential-contact thresholds detailed below |

Every `wt_replicates` entry contains `method`, unique `seed`, `provenance`,
`metrics: {metric_name: value}` and
`binding_measurements: {evaluation_id: {value, metric, unit}}`. The WT envelope
plus preregistered margin is clipped to the preregistered acceptable limit.
WT measurements outside that limit prevent calibration; candidates never relax
it. The declared ligand WT noise bound must cover the observed replicate range.

`spec.binding` contains integer `minimum_replicates >= 2`, integer
`minimum_methods >= 2`, probability `pose_consistency_min`, probability
`contact_reproducibility_min`, and `evaluations` keyed by each evaluation ID.
Every evaluation rule specifies `metric`, `unit`, `wt_variability_bound` and
`binding_support_interval: [low, high]`. Herbicide rules additionally specify
positive `weakening_threshold` and canonical `weakening_direction` (increase
for Kd/Ki/IC50/ΔG; decrease for their supported negative-log quantities).
Native/cofactor evaluations specify `equivalence_interval: [low, high]` for the
WT-relative difference in the same physical quantity/unit.

Every control requires `id`, `kind`, `expected: pass|fail` for structural gate
recovery, `metrics` and `provenance`. The minimum kinds are `resistant`,
`catalytic_loss`, `positive_binding` and `negative_binding`; supplying an empty
allow-list cannot remove them. They also require ligand-resolved measurements:

- Resistant/catalytic-loss controls: `ligand_changes` keyed by evaluation ID,
  each `{value, metric, unit}`. Resistant controls must recover herbicide
  weakening and native equivalence. Catalytic-loss controls must recover native
  impairment.
- Positive/negative binding controls: `ligand_affinities` in the same format.
  Positive controls must fall inside the binding-support interval; negative
  controls must lie beyond its weakened-binding end, respecting endpoint sign.

No control label alone is evidence. Calibration saves source inputs, derived
thresholds, WT variability, control diagnostics, limitations and a digest.
`calibration_errors` recomputes it from saved inputs and rejects modified or
incomplete results. Thresholds must be frozen for a run. Assay feedback can
motivate a **new prospective calibration**, never retroactively approve an old
candidate by editing a bound.

For live mutation generation the same specification also requires `target_id`,
`herbicide_id`, exact `native_ligand_ids`, and `baseline_evidence` entries for
`experimental_structure_comparison`, `published_binding_site`,
`catalytic_residues`, `documented_msa`, `biological_assembly`,
`repeated_structure_prediction` and `repeated_docking`. Each needs
`{validated: true, provenance: ...}`. `neutral_mutation_controls` needs validated
evidence or `{unavailable: true, search_provenance: ...}`. An unavailable
experimental reference is disclosed and keeps this strict live path closed;
the archival diagnostics/report path remains available.

EPSPS additionally requires `evaluated_complexes` containing `glyphosate+S3P`,
`PEP+S3P`, `S3P+PEP`, `apo` and `substrate_bound_control`. Apo is a context
control, not a made-up affinity row.

## Binding assessment record

`assess_binding(record, calibration)` accepts `variant`, `ligand`,
`evaluation_id`, optional `evaluated`, and `methods`. Each method has:

```text
independent_family: string
context: complete mutant/current context
wt_context: corresponding complete wild-type context
replicates: [{seed: integer, artifact: string}, ...]
affinity: affinity object
wt_affinity: affinity object
conclusion: supports_binding | supports_weakening | uncertain
pose_consistency: probability
contact_reproducibility: probability
control_comparison: {
  passed: boolean, positive_control_id: string,
  negative_control_id: string, provenance: scientific provenance
}
```

Independent families must describe genuinely different methods; renamed runs
of the same predictor are not independent methods. Duplicate family/seed
replicates are not counted twice. The importer/run registry owns this identity.
Conclusions are validated against calibrated endpoint intervals and matched WT
comparisons, not simply accepted as annotations.

The seven output states are exactly:

1. `Experimentally confirmed binding`
2. `Experimentally confirmed non-binding or substantially weakened binding`
3. `Strong computational support for binding`
4. `Strong computational support for weakened binding`
5. `Conflicting computational evidence`
6. `Insufficient evidence`
7. `Not evaluated`

No methods and `evaluated: true` is insufficient; no evaluation is not
evaluated. A single pose or a single independent method cannot become strong
support. Equivalent-context method conflict is preserved even when additional
evidence is missing. Different chemistry instead yields non-equivalent evidence.

Optional `experimental` evidence requires experimental provenance,
`context_matched: true`, `direct_binding_assay: true`, an `assay` identifier,
an experimental `measurement` affinity object with Kd or pKd, and `outcome`
`binding` or `weakened_or_nonbinding`. A negative/weakening outcome also requires
`comparison: {wt_affinity, threshold}` with a matched experimental measurement
and a complete interval beyond the threshold. IC50 alone is not a direct
experimental binding assay under this strict contract.

## Candidate screen

`evaluate_candidate(candidate, calibration)` requires `variant`,
`required_evaluations` matching the frozen registry and one `bindings` entry per
evaluation. Each entry contains `evaluation_id`, `ligand`, `role`, `wt_context`,
`mutant_context`, `wt_affinity`, `mutant_affinity` and `assessment` (the input
record above). Its exact WT/mutant affinities and contexts must match one method
record; evidence from another ligand, variant or protocol cannot be reused.

`structural_metrics` maps each of these names to
`{value, unit, provenance}`:

```text
tm_score_forward, tm_score_reverse, ca_lddt, coverage,
global_ca_rmsd, active_site_backbone_rmsd, pocket_sidechain_rmsd,
ligand_rmsd, pocket_volume_change, contact_change, folding_ddg
```

Similarity/coverage use lower bounds; RMSD, stability and perturbations use
upper bounds. Pocket-volume/contact/assembly changes are gated by magnitude.
Assembly-required protocols also require `assembly_change`. Units must match
calibrated units. Scientific applicability exceptions need prospective
protocol design; a missing metric is not an automatic pass.

`biological_context_complete`, `contacts_reproducible` and
`protected_interactions_retained` each require
`{passed: boolean, reason: string, provenance: ...}`.

Output `gates` each have `name`, `status` (`PASS`, `FAIL`, `MISSING` or
`DISAGREEMENT`) and an exact `reason`. Decisions are:

- `METHOD_DISAGREEMENT_REQUIRES_REVIEW` for equivalent-context method conflict.
- `FAILS_COMPUTATIONAL_SCREEN` for a known failed gate (other missing gates remain
  visible).
- `INSUFFICIENT_EVIDENCE` for missing requirements with no established failure.
- `MEETS_COMPUTATIONAL_SCREEN` only when every required gate passes.
- `EXPERIMENTALLY_VALIDATED` / `EXPERIMENTALLY_REJECTED` only for reviewed
  `experimental_validation` records covering scope
  `herbicide_interference_and_native_function`, with experimental provenance,
  `measurement_artifacts` and `outcome: validated|rejected`.

Loading new assay data does not set `reviewed: true`. A computational pass must
never be labelled herbicide resistant.

## Contact and mutation contracts

`extract_pose_contacts(protein_atoms, ligand_atoms, variant=..., ligand=...,
pose_id=..., replicate_id=..., model_id=..., method=...,
canonical_mapping=..., annotations=..., provenance=...)` consumes Cartesian
atom dictionaries (`name`, `element`, `xyz`, and protein `chain`,
`structure_residue`, `amino_acid`). `parse_pdb_atoms(text, model=1)` reads heavy
atoms and preserves insertion codes and the highest-occupancy alternate atom.

General contacts use 4.5 Å; typed hydrogen bonds require 3.5 Å and explicit
D–H–A angle ≥120°; typed salt bridges use 4.0 Å; typed hydrophobic contacts use
4.5 Å; metal coordination hypotheses use 2.8 Å; geometric clash flags use 0.65
of summed van der Waals radii, excluding typed bonded pairs. Untyped coordinates
do not justify hydrogen-bond, aromatic or water-bridge assignments. The
`CONTACT_DEFINITIONS` object exports these definitions and limitations.

`summarize_contacts(rows, pose_registry)` needs **all** sampled poses, including
zero-contact models. Registry entries contain `variant`, `ligand`, `pose_id`,
`replicate_id`, `model_id`, `method`. It exports pose/model frequencies, Wilson
95% intervals over registered independent models, methods, minimum distance,
atoms, types, annotations and provenance. It does not infer true independence.

`classify_residues(residues, native_ligands=[...], thresholds=...)` expects each
residue to contain `canonical_residue`, `protected`, `catalytic`, `conservation`,
`herbicide` and `native: {ligand_id: evidence}`. Each contact evidence object has
`frequency`, `interval`, `method_frequencies: {method: probability}` and integer
`n`. Intervals must contain frequencies. Nomination requires at least two models
and methods, known protection/conservation, all native evidence, reproducibility
and selectivity across the entire interval. Threshold keys are
`reproducibility`, `native_max`, `selectivity_min`, `method_tolerance`,
`conservation_max`, `minimum_models`, `minimum_methods`,
`second_shell_max_angstrom`. Second-shell evidence requires positive
`distance_to_herbicide_contact_angstrom` and scientific `spatial_evidence`, never
sequence adjacency. Independent site support uses `experimentally_reported` or
`cross_method_binding_site_support`.

`validate_mutation_rationale(rationale, classified_residue)` requires `mutation`,
`position_selection`, `herbicide_interaction_disrupted`,
`native_interactions_to_retain`, `structural_support`, `evolutionary_support`,
`expected_physical_effect`, `stability_consequence` with scientific provenance,
`uncertainty`, `falsifying_experiment`, and `conservative_substitution: true`.
The classified residue must also be nomination-eligible. These explanations
organize a falsifiable hypothesis; they do not themselves measure affinity,
stability, retained catalysis or resistance.
