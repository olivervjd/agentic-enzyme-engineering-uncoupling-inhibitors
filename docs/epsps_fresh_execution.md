# Fresh EPSPS execution

The fresh campaign runs EPSPS–glyphosate first. It does not transfer results to the other registered targets. Execution journals distinguish actual tool work, imported evidence, scientific validation and explicit dependency skips.

## Protocol

The prospective input manifest fixes seeds 211, 223 and 227, Boltz 2.2.1, three recycling steps, 100 structure sampling steps, one structure per seed, 200 affinity sampling steps and three internal affinity samples. The latter are not three independent experimental replicates.

The 63 jobs comprise three subjects (mature Arabidopsis WT, E. coli WT and E. coli G96A):

- 36 structure-only jobs: apo, S3P-only, glyphosate with S3P, and PEP with S3P under each seed.
- 27 separate affinity jobs: glyphosate with S3P, PEP with S3P, and S3P scored with PEP under each seed.

Documented ColabFold alignments are reused as verified inputs. Model predictions are generated anew. Arabidopsis residues 77–520 match the observed 7PXY protein. The open apo crystal is compared only to apo predictions. Closed E. coli references retain their homolog scope.

The declared pH 7.4 ligand hypotheses use PEP −3, S3P −3 and glyphosate −2. Dimorphite-DL alternatives are retained. These are modeling assumptions, not measured microstate populations. Boltz's affinity parser changes each scored ligand to a neutral form; its effective chemistry is recorded before inference. Those affinity geometries are excluded from the charged-state contact analysis.

Vina 1.2.7 uses three docking seeds, exhaustiveness 8, at most five poses per seed, and one CPU per worker. Required co-substrates remain rigid receptor components. Docking on predicted receptors and experimental 1G6S redocking remain separate. Scores are empirical Vina scores in kcal/mol, not thermodynamic free energies. Highly charged ligands and protein protonation remain method limitations. Independent seeds, receptor conformations and correlated poses are never pooled into one replicate count.

DDGun3D scores a separately labeled panel of three homolog-mapped diagnostic substitutions and six historical hypotheses on one experimental apo structure. It uses the existing ColabFold alignment through the documented profile/scoring APIs, with narrowly recorded compatibility patches. Raw unfolding ΔΔG and its sign-reversed folding convention are both preserved. No calibrated uncertainty or newly generated mutant structure is claimed.

## Execution

Run module commands from the repository root with the appropriate scientific Python environment. Use new output directories; the runners refuse to overwrite completed campaigns.

```sh
python -m herbicide_desensitization_agent.examples.prepare_epsps_chemical_states --help
python -m herbicide_desensitization_agent.examples.epsps_fresh_campaign prepare --help
python -m herbicide_desensitization_agent.examples.epsps_fresh_campaign run --help
python -m herbicide_desensitization_agent.examples.epsps_fresh_campaign analyze --help
python -m herbicide_desensitization_agent.examples.epsps_vina_parallel --help
python -m herbicide_desensitization_agent.examples.epsps_stability --help
python -m herbicide_desensitization_agent.examples.epsps_evidence_stages --help
python -m herbicide_desensitization_agent.examples.export_fresh_epsps --help
python -m herbicide_desensitization_agent.examples.finalize_fresh_epsps --help
```

Preparation can record a remote MSA root without changing the MSA bytes. Run Boltz on the GPU host, then transfer the effective manifest, original YAMLs, model/affinity/confidence outputs and originating journals. Run postflight analysis on the receiving host to verify the complete cohort, sequence, MSA and output hashes. Do not replace model timings with file modification times.

The parallel docking launcher preserves each subprocess journal and scientific output. Its orchestration timing is separate from scientific-stage timing. The dashboard exporter copies an explicit scientific artifact allowlist, never credential files or a user's home directory.

Evidence synthesis, review and judge use `gpt-5.6-luna`. Local Codex API-key reuse is explicit; credentials are read in memory and are not sent to the compute hosts or written into outputs. Model commentary cannot change numerical results or scientific gates.

## Interpreting completion

A tool can execute successfully and still supply insufficient evidence. The core calibration, nomination-readiness and function gates actually run; missing acceptance bounds or control outcomes cannot be replaced with constants after observing the results. The learning intake check runs, but fitting remains skipped without new traceable candidate assays. No generated dashboard status can override these gates.

The final report distinguishes new calculations, previously published measurements, diagnostic historical substitutions, absent uncertainty and downstream skipped work. A completed computation is not a demonstration of herbicide resistance.
