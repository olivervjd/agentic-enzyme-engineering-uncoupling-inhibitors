# EPSPS evidence and chemical-state protocol

The versioned baseline is `herbicide_desensitization_agent/app/registry/epsps_evidence.json`. It contains curated experimental observations and explicit missing evidence; it does not satisfy the calibration gate by itself.

Arabidopsis EPSPS (AT2G45300/P05466) has a same-target experimental structure, [7PXY](https://www.rcsb.org/structure/7PXY). Its resolved mature sequence, canonical residues 77–520, exactly matches the 444-residue campaign target. It is an open apo structure. Use it to validate matched apo models. A global comparison to a closed complex conflates domain closure with folding error; closed-state assessment needs appropriate homolog structures and domain-aware interpretation. No glyphosate or S3P ligand pose is observed in 7PXY. The crystallization magnesium does not establish a required catalytic metal.

The [Arabidopsis primary study](https://doi.org/10.1016/j.csbj.2022.03.020) measures phosphate-release activity at pH 7.4 in 50 mM HEPES-KOH, 35°C. Standard assays use 1 mM each PEP and S3P; variable-substrate kinetics use 100–1000 µM while the other substrate is fixed at 1 mM. Extraction-buffer EDTA and salt concentrations must not be substituted for final assay concentrations. The pH 7.5 crystal reservoir is a separate condition. The [earlier primary biochemical study](https://doi.org/10.1111/j.1432-1033.1984.tb08378.x) found no evidence of a metal cofactor in a bacterial homolog. Accordingly, the registry retains S3P context and removes the unsupported mandatory EPSPS catalytic-ion requirement. Ion identity, ionic strength and preparation still require auditing.

The [Funke 2009 primary study](https://doi.org/10.1074/jbc.M809771200) supplies quantitative E. coli WT, P101S, T97I and T97I/P101S controls at pH 7.5, 100 mM KCl and 2 mM DTT. Their activity and substrate-utilization tradeoffs remain visible. These are homolog controls, not measured Arabidopsis mutants. Its IC50 measurements use PEP at each enzyme's own Km, so PEP concentration differs across variants. Reported plus-minus errors have no curated confidence-level interpretation. Km, Ki and IC50 are kinetic endpoints and must never be relabeled as Kd. No verified same-target neutral mutant was located in these sources; WT is not a replacement for a neutral mutant.

The canonical resistance positions G177/T178/P182 correspond to E. coli G96/T97/P101; both sequence alignment and the Arabidopsis paper support that mapping. Catalytic K99/D407/E435 map to E. coli K22/D313/E341. Mutagenesis and intermediate-trapping evidence comes from homolog experiments, not experimentally verified plant substitutions. Mature model positions equal canonical positions minus 76; validate identities before mutation.

## Reproducible ligand enumeration

Install the optional `chemistry` dependencies and run:

```bash
python -m herbicide_desensitization_agent.examples.prepare_epsps_chemical_states \
  --input calibration_inputs.json --output chemical_states_ph74.json
```

This runs [Dimorphite-DL](https://doi.org/10.1186/s13321-019-0336-9) 2.0.2 at pH 7.4 with precision 1.0. It retains all enumerated states and verifies covalent identity and stereochemistry with RDKit. The recorded campaign hypotheses are PEP −3, S3P −3 and glyphosate −2 (protonated amine). These are explicit preparation assumptions, not experimentally validated bound microstates or measured populations. No fractions or molecule-specific pKa values are invented. The actual enumeration yielded two PEP, two S3P and four glyphosate alternatives; retain these for charge sensitivity analysis.

Before comparing methods, record the **actual encoded ligand state after each downstream standardization step**. A tool may neutralize a requested charged SMILES; such outputs are not chemically equivalent to docking a charged state. Preserve both requested and encoded identities, net charges, software version and transformation. Do not call matching request strings evidence of matching chemistry. Keep the same prepared state across WT and variants. Revisit conclusions that depend on an untested microstate; do not select a best-scoring state after inspecting candidate outcomes.

Protein protonation, local pKa shifts, catalytic waters and assay ionic strength remain separate questions. Ligand-only enumeration does not complete their validation. Calibration must still recover resistance and native-function tradeoffs with appropriate controls and uncertainty before nomination.
