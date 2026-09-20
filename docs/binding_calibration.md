# Binding Calibration

The implementation is operational, but the current EPSPS dataset has **no
endpoint-matched experimental labels**. Its report returns `INSUFFICIENT_DATA`,
leaves the fitted model null, and keeps function-retention gates closed.
No measurements or calibration coefficients are invented for the demo.

## Method and Limits

- Affine fit: Boltz predicted pIC50 to measured pIC50, within one target, assay,
  chemical-series and prediction-protocol scope.
- Exact IC50 measurements in M, mM, uM or nM are converted to -log10 molar. Kd,
  Ki, Km, censored values and mismatched scopes are excluded.
- Predeclared independent splits: at least 8 training, 9 uncertainty-calibration,
  and 10 test groups. Aggregate seed/assay replicates before splitting. A subject
  cannot cross splits even with different ligands.
- Held-out absolute residuals define a 90% split-conformal interval using the
  `ceil((n + 1) * 0.9)` order statistic. Test observations never fit coefficients.
- Test gates: positive slope; RMSE <= 0.5 pIC50 and no worse than raw predictions;
  coverage >= 90%; interval half-width <= 1 pIC50; held-out inputs in training range.
- No extrapolation or application to a different scope. Source/prediction hashes
  are verified and the application ensemble is bound by hash.

These minimum counts and limits are **engineering safeguards, not a statistical
power calculation**. Nine calibration groups merely permit a finite 90% interval.
Correlated mutation families require appropriate grouping. The standard
[split-conformal framework](https://arxiv.org/abs/2210.14735) assumes exchangeability;
coverage is marginal, not guaranteed for arbitrary new mutations or every subgroup.
Do not repeatedly tune against the test set. Passing IC50 calibration does not
establish Kd, native-substrate equivalence, catalytic turnover or plant resistance.

## Real Control Outcome

Matched E. coli G96A seeds give mean mutant-minus-WT predicted pIC50 differences
of about -0.034 for glyphosate and +0.232 for PEP. The glyphosate direction agrees
qualitatively with insensitivity; PEP conflicts with the reported weaker affinity.
These are warnings, not numeric calibration labels or Arabidopsis validation.

Sources: [1MI4 and primary study](https://www.rcsb.org/structure/1MI4),
[official Boltz endpoint definition](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).

## Run

```sh
.venv/bin/python -m herbicide_desensitization_agent.examples.calibrate_binding \
  --calibration /path/to/recorded-calibration \
  --output /path/to/binding_calibration.json
```

Exit code 2 means evaluated but not validated. Add `--dataset` for measured data.
The integrated workflow accepts `--binding-dataset`, writes its report before the
structure-quality gate, and cannot use IC50 calibration to override function gates.
New structure-calibration runs generate a binding report automatically.

## Dataset Contract

Top-level fields: `application_prediction_sha256` (current `predictions.json` hash),
`scope`, and `observations`. The nonempty scope fields are `target_id`, `assay_id`,
`conditions` (pH, temperature, substrate concentrations and other assay details),
`prediction_protocol` (pinned model, chemistry, MSA and settings), and `ligand_scope`.

Every observation requires:

- Unique `id`, `subject`, `ligand`, independent `group`, `split` (`train`,
  `calibration`, or `test`), and the exact same `scope` object.
- `endpoint: "IC50"`, `relation: "="`, positive numeric `value`, `unit`, and
  `evidence_type: "experimental"`.
- `source`: HTTPS `url`, table/page `locator`, relative local `path` and `sha256`.
- `prediction`: relative local `path`, `sha256`, and unique `record_ids`. The file
  must contain `records` with matching `id`, `subject`, `ligand`, and numeric
  `predicted_pIC50`. The loader derives the mean rather than trusting a supplied mean.

Keep artifacts inside the dataset directory. Hashes prove byte identity, not
scientific truth: transcription, mutation numbering, assay comparability and
protocol identity require curation. Every report includes a scientific legend.
Synthetic software tests are never exported as scientific demo evidence.
