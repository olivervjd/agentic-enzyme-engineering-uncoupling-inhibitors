# Hackathon Demo

A local, read-only workbench with actual recorded results. It does not simulate
a live pipeline, launch GPU jobs or approve candidates. Runtime libraries and data
are bundled locally; the built demo makes no model API requests.

## Launch

From the repository root:

```sh
.venv/bin/python -m herbicide_desensitization_agent.examples.prepare_demo \
  --calibration /path/to/epsps-calibration-matchedmsa-20260920 \
  --archive /path/to/full-workflow-20260919-docking-affinity \
  --model-access /path/to/model_access.json \
  --output demo/public/data
cd demo
npm ci
npm run dev
```

`--model-access` is optional. Export verifies calibration artifact hashes, parses
coordinates with Biopython and keeps older candidates in a separately labeled
archive. Generated data are Git-ignored. The server binds to loopback and serves
only the demo root. For a local production bundle: `npm run build`, then
`npm run preview`. Do not place secret files in the public directory.

## Demo Sequence

1. Show EPSPS, focus the ligand pocket, change ligand context and seed, and overlay
   another seed of the same protein/context, aligned by C-alpha least squares.
2. Contrast 30 passing structural comparisons with missing binding calibration.
   Stable predicted shape is not proof of function.
3. Open Binding Calibration and explain the G96A PEP direction warning.
4. Open Candidate Archive: six earlier mutations plus native WT, with the older
   query-only MSA protocol and high WT variability clearly disclosed.
5. Open Evidence & Review: sources, a real Luna calibration review, and blocked
   workflow stages. GPT-Rosalind was not listed for the current key.

Tables/figures have expandable legends. Other registered targets show no-run states,
not fabricated results. A standalone model review is not a completed candidate judge.

## Codex API Key

Explicit opt-in reuses only an actual stored API key in Python memory:

```sh
.venv/bin/python -m herbicide_desensitization_agent.examples.check_model_access \
  --use-codex-api-key --judge-model gpt-5.6-luna \
  --binding-report /path/to/binding_calibration.json \
  --output /path/to/model_access.json
```

The integrated workflow accepts the same key flag and separate evidence/review/judge
model IDs. `OPENAI_API_KEY` takes precedence. ChatGPT OAuth tokens are never used as
API credentials. No automatic model substitution occurs. Catalog listing and a
verified inference request are reported separately. See
[OpenAI authentication](https://learn.chatgpt.com/docs/auth).

## Verification

```sh
.venv/bin/python -m unittest discover -s tests
cd demo
npm run build
npx playwright install chromium
npm run test:browser
```

Use `PLAYWRIGHT_CHANNEL=chrome` for installed Chrome, `DEMO_URL` for a different
port, and `DEMO_SCREENSHOTS` for an output directory. Tests use fresh contexts,
check canvas pixels and interactions, and capture desktop/mobile screenshots.
When the sandbox blocks launching Chrome, use the in-app browser's Playwright and
screenshot APIs instead. Do not change security settings or use personal profiles.
