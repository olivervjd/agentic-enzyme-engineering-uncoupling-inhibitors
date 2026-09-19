# Herbicide Desensitization Agent

Milestone 1 prototype for a review-gated computational decision engine covering
six fixed *Arabidopsis thaliana* AGI–herbicide pairings across five herbicide
classes.

The current implementation contains a fixed provenance-aware target registry,
typed schemas, deterministic validators, mock scientific backends, an end-to-end
orchestrator, six placeholder examples, and tests.

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
```

Real NVIDIA BioNeMo/NIM and GPT-Rosalind clients should be added behind the
interfaces in `app/backends/`; the orchestrator does not depend on vendor SDKs.
