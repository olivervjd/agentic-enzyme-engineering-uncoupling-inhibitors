# G177 saturation scan and Foldseek comparison

These examples implement the explicitly requested 19 substitutions at Arabidopsis
EPSPS P05466 G177, with wild type as the twentieth sequence. They do not change
automatic nomination rules or approve a mutation. Canonical residues 77–520 form
the 444-residue mature construct; canonical G177 is model residue 101.

Two separate campaigns are supported:

- `g177_scan`: PEP + S3P and glyphosate + S3P cofolding with Boltz2 affinity
  scoring, for three seeds per condition (120 predictions).
- `g177_apo_foldseek`: fresh protein-only predictions for the same 20 sequences
  and three seeds (60 predictions), followed by Foldseek comparisons against
  all three newly predicted WT structures (180 comparisons).

The apo table summarizes the three comparisons with matching seeds. WT
self-comparisons are shown separately from cross-seed WT variability. The
reference is a prediction, not an experimental structure.

## Environment and inputs

Use a GPU environment with standalone `boltz==2.2.1`, its downloaded model cache,
Biopython and RDKit. The recorded runs used Python 3.12 on a Brev NVIDIA H200.
Boltz is run through its CLI; these commands do not use NVIDIA NIM endpoints.
Install this project with its `calibration` and `reports` extras where needed:

```bash
python -m pip install -e '.[calibration,reports]'
```

Install [Foldseek](https://github.com/steineggerlab/foldseek#installation) on the
machine performing structure comparisons. No Foldseek sequence or structure
database download is needed for these explicit pairwise comparisons.

The affinity preparation step requires a completed
[fresh EPSPS campaign](epsps_fresh_execution.md), containing
`calibration_inputs.json`, its `arabidopsis_wt` subject and registered A3M file.
It verifies the canonical sequence and donor MSA hashes. The same homolog rows
are retained across mutants; only the query residue changes. Model weights,
MSAs, generated structures and local results are runtime inputs/outputs, not
bundled with this code.

Set paths for your existing environment:

```bash
G177_DONOR=/path/to/epsps-fresh-campaign
G177_AFFINITY=/path/to/new-g177-affinity
G177_APO=/path/to/new-g177-apo
G177_BOLTZ=/path/to/boltz-venv/bin/boltz
G177_CACHE=/path/to/boltz-cache
G177_FOLDSEEK=/path/to/foldseek/bin/foldseek
G177_REPORT=/path/to/g177-structure-report
```

Preparation requires a fresh campaign directory. Execution refuses duplicate
launches. Run commands in the Boltz environment, without Python's `-O` option:
integrity assertions must remain enabled.

## Ligand-bound cofolding and affinity

```bash
python -m herbicide_desensitization_agent.examples.g177_scan prepare \
  --root "$G177_AFFINITY" --donor "$G177_DONOR"
python -m herbicide_desensitization_agent.examples.g177_scan run \
  --root "$G177_AFFINITY" --cache "$G177_CACHE" --boltz "$G177_BOLTZ"
python -m herbicide_desensitization_agent.examples.verify_g177_scan \
  "$G177_AFFINITY" > "$G177_AFFINITY/validation.json"
python -m herbicide_desensitization_agent.examples.plot_g177 \
  --results "$G177_AFFINITY/results.json" --output-dir "$G177_AFFINITY/report"
python -m herbicide_desensitization_agent.examples.package_g177_scan \
  "$G177_AFFINITY" "${G177_AFFINITY}-artifacts.tar.gz"
```

`collect --root ...` refreshes partial or final affinity results without launching
inference. Outputs include the sequence FASTA, input/hash manifest, execution
logs, raw predictions, `results.json` and `binding_table.md`. The plot command
exports PNG and SVG. Packaging requires a complete, audited campaign.

Affinity results retain the raw `affinity_pred_value` scale, log10(IC50 / µM),
where lower scores predict stronger binding. These are not measured Kd or PEP
Km. The Boltz affinity parser neutralizes the scored ligand; the manifest records
both requested and encoded chemistry. S3P is present in both ligand contexts.

## Fresh apo structures and Foldseek

The apo preparation uses the affinity campaign's sequence/MSA manifest; it does
not reuse ligand-bound coordinates. Preparing the affinity inputs is sufficient
if only the apo comparison is wanted.

```bash
python -m herbicide_desensitization_agent.examples.g177_apo_foldseek prepare \
  --root "$G177_APO" --donor "$G177_AFFINITY"
python -m herbicide_desensitization_agent.examples.g177_apo_foldseek run \
  --root "$G177_APO" --cache "$G177_CACHE" --boltz "$G177_BOLTZ"
python -m herbicide_desensitization_agent.examples.g177_apo_foldseek analyze \
  --root "$G177_APO" --foldseek "$G177_FOLDSEEK"
python -m herbicide_desensitization_agent.examples.g177_apo_foldseek pack \
  --root "$G177_APO"
python -m herbicide_desensitization_agent.examples.report_g177_structure \
  --results "$G177_APO/results.json" --output-dir "$G177_REPORT"
cp "${G177_APO}-artifacts.tar.gz" "$G177_REPORT/structures-and-foldseek.tar.gz"
```

The report contains a Markdown table, PNG/SVG column plots and a copy of the full
results JSON. The archive contains original mmCIFs, protein-only PDB exports,
confidence files, requests, manifests, logs and raw Foldseek alignments. Apo MSA
contents are omitted from the archive; their hashes and paths remain recorded.

Both campaigns use seeds 211, 223 and 227, three recycles, 100 structure sampling
steps and one structure sample per seed. The affinity campaign additionally uses
200 affinity sampling steps and three internal affinity samples. Results are
reported as means and sample SDs, not calibrated confidence intervals.

Foldseek runs `easy-search --alignment-type 1 --exhaustive-search 1 -a 1` to
obtain WT-length-normalized TM-score, structural lDDT and alignment spans.
Structural lDDT is a coordinate comparison, not prediction-confidence pLDDT.
`alnlen` includes gap columns; `paired_residues` explicitly excludes them, and
the table reports the range of paired positions across seeds. Biopython computes
a separate least-squares Cα RMSD using all 444 sequence-matched positions.

Every structure is checked for the exact sequence, residue 101, protein-only
chain identity and complete Cα coordinates. The run checks unique job counts,
hashes, metric ranges and raw alignment strings. Packing rechecks coordinate
hashes to reject artifacts changed since analysis.

## Recorded validation and interpretation

On 20 September 2026, all 120 ligand-bound and 60 fresh apo predictions completed
on the H200. All 180 apo/WT Foldseek comparisons were obtained using Foldseek
`463739e0014a1549a527de589102cde98f802f37`. G177S had the closest mean matched-seed
apo structure (TM-score 0.9973; full-chain Cα RMSD 0.422 Å). WT predictions varied
by 0.344–3.524 Å across distinct seeds. That variability limits attribution of
predicted differences to mutation alone; none of these scores establishes
retained PEP binding, glyphosate resistance, stability or catalytic activity.

The sequence/MSA and gapped-alignment regression tests use synthetic fixtures:

```bash
python -m unittest discover -s tests -p 'test_g177_scan.py' -v
```
