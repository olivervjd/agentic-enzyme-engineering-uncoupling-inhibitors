"""Render the complete validated G177 apo structure comparison."""
import argparse
import copy
import json
from pathlib import Path
import sys
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

if sys.version_info[:2] == (3,14) and matplotlib.__version__ == '3.10.1':
    from matplotlib.path import Path as MplPath
    def copy_path(self, memo=None):
        result = copy.copy(self)
        if memo is not None: memo[id(self)] = result
        result._vertices = self.vertices.copy()
        result._codes = None if self.codes is None else self.codes.copy()
        result._readonly = False
        return result
    MplPath.__deepcopy__ = copy_path
    MplPath.deepcopy = copy_path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--results', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
source_text = args.results.read_text()
data = json.loads(source_text)
assert data['complete'] and data['predictions']==60 and data['foldseek_comparisons']==180
rows = data['rows']
assert len(rows)==20 and rows[0]['variant']=='WT'
(OUT/'results.json').write_text(source_text)
for row in rows:
    for key in ['ttmscore','lddt','qcov','tcov','full_chain_ca_rmsd_angstrom']:
        assert row[key]['n']==3 and all(math.isfinite(row[key][s]) for s in ['mean','sd','min','max'])

def score(row, key, digits=3):
    x = row[key]
    return f"{x['mean']:.{digits}f} ± {x['sd']:.{digits}f}"

baseline = data['wt_cross_seed_summary']
lines = ['# EPSPS G177: fresh apo refolding and Foldseek comparison', '',
    'All 19 substitutions and wild type were predicted afresh with standalone Boltz 2.2.1. Three seeds (211, 223, 227) were used per sequence: 60 protein-only structures total. No ligands, structure templates, or affinity predictions were used in this run.', '',
    'The target is Arabidopsis EPSPS P05466, mature residues 77–520 (444 residues). Canonical position 177 corresponds to model position 101. Each mutant was compared with the newly predicted wild type from the same seed. The reference is predicted, not experimentally determined.', '',
    '| Amino acid at 177 | Variant | Foldseek TM-score ↑ | Foldseek structural lDDT ↑ | Full-chain Cα RMSD (Å) ↓ | Paired residues (range / 444) |',
    '|---|---|---:|---:|---:|---:|']
for row in rows:
    lines.append(f"| {row['amino_acid']} — {row['name']} | {row['variant']} | {score(row,'ttmscore')} | {score(row,'lddt')} | {score(row,'full_chain_ca_rmsd_angstrom')} | {row['paired_residues']['min']}–{row['paired_residues']['max']} / 444 |")
lines += ['', 'Values are mean ± sample SD across three matched-seed comparisons. Higher TM-score and structural lDDT indicate greater structural similarity; lower RMSD indicates closer Cα positions after superposition. TM-score is normalized by WT length. Structural lDDT is a comparison between coordinates, not the model’s pLDDT confidence.', '',
    '**Wild-type controls.** The WT row compares each WT structure to itself and therefore has perfect scores by construction. Across the three distinct WT seed pairs, '+
    f"TM-score was {score(baseline,'ttmscore')}, structural lDDT {score(baseline,'lddt')}, and full-chain Cα RMSD {score(baseline,'full_chain_ca_rmsd_angstrom')} Å "+
    f"(range {baseline['full_chain_ca_rmsd_angstrom']['min']:.3f}–{baseline['full_chain_ca_rmsd_angstrom']['max']:.3f} Å). These three comparisons share structures and are descriptive controls, not independent experimental replicates.", '',
    '## Method', '',
    '- Standalone Boltz 2.2.1: 3 recycling steps, 100 sampling steps, 1 diffusion sample per seed. Identical homolog MSA rows were reused across variants, changing only the query residue.',
    f"- Foldseek version `{data['foldseek_version']}`: `easy-search --alignment-type 1 --exhaustive-search 1 -a 1`. This uses TM-align global structural alignment. All 60 structures were compared with all three WT structures (180 comparisons); the table uses the 60 matched-seed comparisons.",
    '- Protein-only PDB coordinates were exported from the original mmCIFs. Foldseek supplies TM-score, structural lDDT and alignment coverage. Biopython separately computes least-squares Cα RMSD over all 444 sequence-matched residues, regardless of Foldseek alignment trimming.',
    '- Validation checked every output sequence, residue 177, chain identity, all 444 Cα atoms, input/MSA hashes, comparison uniqueness and finite metric ranges. Raw alignment strings were checked separately: paired residues count only positions with residues in both structures, excluding gap columns.', '',
    '## Interpretation limits', '', *['- '+s for s in data['limitations']], '',
    'A high global similarity score does not rule out a local binding-pocket change. These comparisons cannot establish retained PEP binding, reduced glyphosate binding, protein stability, or catalytic activity.', '',
    '## Files and sources', '',
    '- `results.json`: aggregate metrics, every comparison, WT controls, commands and artifact hashes.',
    '- `structures-and-foldseek.tar.gz`: original mmCIFs, protein PDBs, confidence files, requests, manifests and raw Foldseek output. MSA contents are omitted from the archive; their source and hashes are recorded.',
    '- [Foldseek documentation](https://github.com/steineggerlab/foldseek#alignment-mode)',
    '- [Boltz prediction documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)', '']
(OUT/'comparison-table.md').write_text('\n'.join(lines))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
    'axes.spines.right':False,'axes.spines.left':False,'axes.spines.bottom':False,
    'svg.fonttype':'none','savefig.facecolor':'white'})
fig, axes = plt.subplots(2,1,figsize=(15,10),sharex=True)
fig.subplots_adjust(left=.080,right=.985,bottom=.19,top=.835,hspace=.30)
x = np.arange(len(rows))
colors = ['#97A7B3'] + ['#2878AD']*19
error = dict(ecolor='#263A49',elinewidth=.8,capsize=2,capthick=.8)
for ax in axes:
    ax.axvspan(-.48,.48,color='#E9EEF2',zorder=0)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True,color='#E3E8EC',linewidth=.75)
    ax.tick_params(axis='both',length=0,pad=8,colors='#344554')
    ax.set_xlim(-.7,len(rows)-.3)
width=.36
for offset,key,color,label in [(-width/2,'ttmscore','#2878AD','TM-score'),(width/2,'lddt','#DC8A35','Structural lDDT')]:
    axes[0].bar(x+offset,[r[key]['mean'] for r in rows],width=width,
        yerr=[r[key]['sd'] for r in rows],color=color,label=label,error_kw=error,zorder=3)
axes[0].set_ylim(0,1.08)
axes[0].set_yticks(np.arange(0,1.01,.2))
axes[0].set_ylabel('Structural similarity\nHigher = closer to wild type',labelpad=12)
axes[0].legend(loc='upper right',bbox_to_anchor=(1,1.15),ncol=2,frameon=False)
axes[1].bar(x,[r['full_chain_ca_rmsd_angstrom']['mean'] for r in rows],width=.68,
    yerr=[r['full_chain_ca_rmsd_angstrom']['sd'] for r in rows],color=colors,error_kw=error,zorder=3)
wt = baseline['full_chain_ca_rmsd_angstrom']
axes[1].axhspan(wt['min'],wt['max'],color='#D8B361',alpha=.28,zorder=1,
    label=f"WT across seeds: {wt['min']:.2f}–{wt['max']:.2f} Å")
axes[1].set_ylim(bottom=0)
axes[1].set_ylabel('Full-chain Cα RMSD (Å)\nLower = closer to wild type',labelpad=12)
axes[1].legend(loc='upper right',bbox_to_anchor=(1,1.14),frameon=False)
axes[1].set_xticks(x,['WT (G)' if r['variant']=='WT' else r['variant'] for r in rows],
    rotation=45,ha='right',rotation_mode='anchor')
axes[1].set_xlabel('Substitution at canonical residue 177',labelpad=10)
fig.text(.080,.948,'EPSPS G177: structure after refolding',fontsize=22,fontweight='bold',color='#203444')
fig.text(.080,.908,'Fresh Boltz2 protein-only predictions compared with predicted wild type using Foldseek',fontsize=12,color='#465C6D')
fig.text(.080,.090,'Bars: mean ± sample SD of 3 matched-seed comparisons. WT bars are self-comparisons; the gold band shows WT sampling variation.',fontsize=9.5,color='#465C6D')
fig.text(.080,.057,'Mature EPSPS: 444 residues. Structural similarity does not establish retained binding, stability or enzyme activity.',fontsize=9.5,color='#465C6D')
for ext in ['png','svg']:
    fig.savefig(OUT/f'column-plot.{ext}',dpi=200)
plt.close(fig)
print(json.dumps({'rows':len(rows),'wt_control_rmsd':wt,
    'mutant_tm_range':[min(r['ttmscore']['mean'] for r in rows[1:]), max(r['ttmscore']['mean'] for r in rows[1:])],
    'mutant_rmsd_range':[min(r['full_chain_ca_rmsd_angstrom']['mean'] for r in rows[1:]),max(r['full_chain_ca_rmsd_angstrom']['mean'] for r in rows[1:])]},indent=2))
