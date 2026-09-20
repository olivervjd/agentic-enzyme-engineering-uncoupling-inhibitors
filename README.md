# Agentic Enzyme Engineering: Uncoupling Inhibitors

A research prototype for studying protein substitutions that could reduce
herbicide binding while preserving native enzyme function. It combines an agent
workflow for evidence, analysis and review with reproducible scientific tools,
recorded inputs, structural comparisons and reports.

The main implemented case study is **Arabidopsis EPSPS and glyphosate**. The
workflow compares the native substrate phosphoenolpyruvate (PEP) with glyphosate,
retaining shikimate-3-phosphate (S3P) where the modeled binding context requires
it. The latest examples scan all 19 substitutions at **G177**, predict structures
and affinity scores with **Boltz2**, and compare freshly predicted ligand-free
mutants with wild type using **Foldseek**.

This repository was previously named `herbicide-resistance-workflow`. The Python
package is still named `herbicide_desensitization_agent`.

## What is included

- An orchestrator connecting evidence retrieval, structure-quality checks,
  ligand-contact analysis, mutation nomination, scoring, review and assay intake.
- Control and calibration checks that block automatic nomination when evidence
  is missing or insufficient. Explicit G177 scan commands are separate from
  automatic nomination and do not change its rules.
- Boltz2 structure/affinity campaigns, Foldseek comparisons and optional
  docking/stability diagnostics.
- JSON manifests, hashes, execution records, Markdown tables, PNG/SVG plots and
  archives of scientific artifacts.
- A local browser workbench for inspecting previously generated evidence.

The registry contains six Arabidopsis protein–herbicide pairings: EPSPS–glyphosate,
GLN1;2–glufosinate, GLN2–glufosinate, PsbA–atrazine, ALS/AHAS–chlorsulfuron and
TIR1–2,4-D. Recorded EPSPS results do not validate the other targets.

**Research status:** model scores and structural similarity are computational
evidence, not proof of herbicide resistance or retained catalytic activity.
Synthetic examples are explicitly labeled. Experimental validation remains a
separate requirement.

## Requirements

| Use | Requirements |
|---|---|
| Synthetic example | Git, Python 3.11+ and pip; no GPU or API key |
| Scientific analysis, plotting and full Python test suite | Python 3.12 is the tested environment; install the optional dependencies below |
| Fresh Boltz2 predictions | A separate Linux GPU environment with a compatible NVIDIA driver/PyTorch CUDA setup, `boltz==2.2.1`, model weights/cache and prepared input MSAs |
| Foldseek comparisons | The Foldseek executable; this small comparison set can run on CPU without a downloaded structure database |
| Browser workbench | Node.js and npm; use Node 24.12+ or 22.20+ for the locked dependencies |
| Live explanatory model roles | `OPENAI_API_KEY` and access to the configured evidence, reviewer and judge models |
| Optional advanced stages | Separately configured DiffDock, AutoDock Vina or DDGun3D, depending on the selected workflow |

The recorded G177 GPU runs used an NVIDIA H200 through Brev. An H200 is the
tested machine, not a declared minimum hardware requirement. Model weights,
MSAs and completed campaign artifacts are not installed by this package.

## Installation

The commands below use a macOS/Linux shell. Run Python module commands from the
repository root.

```bash
git clone https://github.com/oliver-intransition/agentic-enzyme-engineering-uncoupling-inhibitors.git
cd agentic-enzyme-engineering-uncoupling-inhibitors

python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[structure,calibration,chemistry,reports]'
```

The optional dependency groups are defined in [pyproject.toml](pyproject.toml):

| Extra | Packages / purpose |
|---|---|
| `structure` | Biopython, NumPy and TM-tools for coordinate analysis |
| `calibration` | Structure dependencies plus RDKit for calibration and chemical checks |
| `chemistry` | Dimorphite-DL and RDKit for ligand-state preparation |
| `reports` | Matplotlib for figures |

For only the synthetic example, `python -m pip install -e .` is sufficient.
The complete test suite needs the scientific extras shown above.

### GPU environment

Set up Boltz separately from any BioNeMo environment. On the GPU host, from this
repository's root, a dedicated environment can be created with:

```bash
python3.12 -m venv work/boltz-venv
source work/boltz-venv/bin/activate
python -m pip install --upgrade pip
python -m pip install 'boltz==2.2.1' -e '.[calibration,chemistry,reports]'
boltz predict --help
```

Ensure CUDA is available to PyTorch and provision the Boltz model cache before
starting a campaign. Follow the [Boltz installation and prediction documentation](https://github.com/jwohlwend/boltz)
for the host-specific setup. Install [Foldseek](https://github.com/steineggerlab/foldseek#installation)
on the machine doing structural comparisons and verify it with `foldseek version`.
The G177 examples use standalone Boltz, not a NVIDIA NIM endpoint.

## Run a quick example

```bash
python -m herbicide_desensitization_agent.examples.run_all
```

This runs the six registry pairings through **synthetic backends** and prints a
JSON summary with `synthetic: true`. It checks the orchestration without network
access, GPU inference or real binding predictions.

## Run the G177 scan

Use the [G177 saturation-scan guide](docs/g177_saturation_scan.md) for the complete
commands, input format, validation and packaging steps.

The scan starts from an existing EPSPS input campaign containing
`calibration_inputs.json`, its registered Arabidopsis WT sequence and a verified
ColabFold A3M alignment. A fresh clone alone does not contain these inputs. The
preparation step verifies sequence/MSA hashes and changes only the query residue,
keeping the homolog alignment rows constant.

There are two independent prediction stages:

| Stage | Predictions | Main results |
|---|---:|---|
| Ligand-bound scan | 20 sequences × 2 ligand contexts × 3 seeds = 120 | PEP/glyphosate affinity scores, structures and comparison plots |
| Fresh apo refolding | 20 sequences × 3 seeds = 60 | Protein-only structures; 180 Foldseek comparisons against three predicted WT references |

For example, after preparing an affinity input campaign, run fresh apo refolding
and generate its report in the Boltz environment:

```bash
python -m herbicide_desensitization_agent.examples.g177_apo_foldseek prepare \
  --root /path/to/new-g177-apo \
  --donor /path/to/g177-affinity-campaign

python -m herbicide_desensitization_agent.examples.g177_apo_foldseek run \
  --root /path/to/new-g177-apo \
  --cache /path/to/boltz-cache \
  --boltz /absolute/path/to/boltz-venv/bin/boltz

python -m herbicide_desensitization_agent.examples.g177_apo_foldseek analyze \
  --root /path/to/new-g177-apo \
  --foldseek /absolute/path/to/foldseek

python -m herbicide_desensitization_agent.examples.report_g177_structure \
  --results /path/to/new-g177-apo/results.json \
  --output-dir outputs/g177-structure-report
```

Replace every `/path/to/...` value with a real path. Use a new campaign directory;
preparation and inference reject duplicate launches. Do not use Python's `-O`
flag, because the campaign integrity assertions must remain enabled.

The report summarizes TM-score, structural lDDT, paired-residue counts and
full-chain Cα RMSD, with mean ± SD across three matching seeds. It also reports
WT variation across seeds. Canonical residue 177 is residue 101 in the mature
444-residue construct. The WT comparator in this scan is predicted, not experimental.

## Run the integrated EPSPS agent workflow

The integrated workflow requires an existing precursor campaign, a Boltz runtime
and its cache. See the [live EPSPS guide](docs/live_epsps.md) for those input
contracts, reuse options and review behavior, and the
[fresh execution guide](docs/epsps_fresh_execution.md) for the expanded scientific
campaign stages.

```bash
python -m herbicide_desensitization_agent.examples.run_integrated_epsps \
  --previous /path/to/completed-docking-affinity-campaign \
  --calibration /path/to/new-calibration \
  --output /path/to/new-agent-workflow \
  --boltz-cache /path/to/boltz-cache \
  --boltz /absolute/path/to/boltz-venv/bin/boltz \
  --evidence-calibration /path/to/reviewed-frozen-calibration.json
```

Supply `OPENAI_API_KEY` through the process environment for live explanatory
model calls. `EVIDENCE_MODEL`, `REVIEW_MODEL` and `JUDGE_MODEL`, or their CLI
options, select the model IDs. Use IDs available to your account. The G177
scientific scan and deterministic reports do not require an OpenAI key.

Useful integrated-workflow options:

- `--diagnostics-only`: inspect an existing calibration without Boltz inference
  or mutation generation; model review/retrieval can still make network calls.
- `--cached-only`: require existing predictions; scientific gates still apply.
- `--curated-evidence-only`: disable model-based evidence synthesis explicitly.
- `--no-literature-retrieval`: disable online Europe PMC retrieval.
- `--assays /path/to/assays.json`: supply real, traceable assay results for learning.

Live mutation nomination requires a reviewed, prospectively frozen calibration.
Missing evidence blocks nomination. Exit code **2** means `NEEDS_REASSESSMENT`;
inspect `execution.json` and `workflow_stages.json` for the reason.

## Run the browser workbench

The workbench displays recorded results; it does not launch GPU jobs. First
export your completed campaign data:

```bash
python -m herbicide_desensitization_agent.examples.prepare_demo \
  --calibration /path/to/completed-calibration \
  --archive /path/to/recorded-mutation-campaign \
  --workflow-root /path/to/completed-agent-workflow \
  --output demo/public/data

cd demo
npm ci
npm run dev
```

Open the local URL printed by Vite. Use `npm run build` and `npm run preview`
for a production build. Generated workbench data are Git-ignored and must be
exported before the UI can display a campaign. See the
[workbench guide](docs/evidence_workbench.md) for a self-contained offline report.

## Tests

From the repository root, with the scientific extras installed:

```bash
python -m unittest discover -s tests -v
python -m unittest discover -s tests -p 'test_g177_scan.py' -v
```

The second command runs just the G177 sequence, MSA, alignment and archive
integrity tests. Test fixtures do not establish biological accuracy. Browser
test setup is described in the [demo guide](docs/hackathon_demo.md#verification).

## Repository layout and further reading

| Path | Contents |
|---|---|
| `herbicide_desensitization_agent/app/` | Agents, orchestration, backends, schemas, validators and target registry |
| `herbicide_desensitization_agent/examples/` | Executable campaign, analysis, export and reporting modules |
| `tests/` | Python regression and integrity tests |
| `demo/` | Vite-based molecular evidence workbench |
| `docs/` | Scientific protocols, input contracts and detailed run instructions |
| `outputs/` | Local generated artifacts; Git-ignored |

- [G177 saturation scan and Foldseek comparison](docs/g177_saturation_scan.md)
- [Control-first live EPSPS workflow](docs/live_epsps.md)
- [Fresh EPSPS execution](docs/epsps_fresh_execution.md)
- [Ligand chemical-state protocol](docs/epsps_chemical_protocol.md)
- [Binding calibration](docs/binding_calibration.md)
- [Evidence contracts](docs/evidence_contracts.md)
- [Evidence workbench](docs/evidence_workbench.md)

Affinity-model scores, docking scores, structural confidence and experimental
binding/function measurements are different quantities. In particular, a
similar predicted fold does not demonstrate preserved PEP binding or reduced
glyphosate binding. Keep the recorded endpoint, chemical context and uncertainty
with every reported result.
