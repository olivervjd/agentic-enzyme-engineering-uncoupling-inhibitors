"""Plot a completed G177 ligand-affinity scan without rerunning inference."""
import argparse
import json
import copy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Matplotlib 3.10.1's deepcopy(super()) recurses under Python 3.14.
# Apply a process-local equivalent that copies Path arrays without that call.
if sys.version_info[:2] == (3, 14) and matplotlib.__version__ == '3.10.1':
    from matplotlib.path import Path as MplPath
    def copy_path(self, memo=None):
        result = copy.copy(self)
        if memo is not None:
            memo[id(self)] = result
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
args.output_dir.mkdir(parents=True, exist_ok=True)
data = json.loads(args.results.read_text())
assert data['complete'] and data['completed'] == 120
rows = data['rows']
assert len(rows) == 20 and all(r[k]['n'] == 3 for r in rows for k in ('pep', 'glyphosate'))

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.spines.left': False, 'axes.spines.bottom': False,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white',
})
fig, ax = plt.subplots(figsize=(15, 7.4), dpi=150)
fig.subplots_adjust(left=.078, right=.985, bottom=.245, top=.77)
x = np.arange(len(rows))
width = .36
ax.axvspan(-.48, .48, color='#E9EEF2', zorder=0)
for offset, key, color, label in (
    (-width/2, 'pep', '#2878AD', 'PEP'),
    (width/2, 'glyphosate', '#DC8A35', 'Glyphosate'),
):
    means = [row[key]['mean'] for row in rows]
    sd = [row[key]['sd'] for row in rows]
    ax.bar(x + offset, means, width=width, yerr=sd, color=color,
           linewidth=0, label=label, zorder=3,
           error_kw={'ecolor': '#303B44', 'elinewidth': .9,
                     'capsize': 2.4, 'capthick': .9})

ax.set_xticks(x, ['WT (G)' if row['variant'] == 'WT' else row['variant'] for row in rows],
              rotation=45, ha='right', rotation_mode='anchor')
ax.set_xlim(-.7, len(rows)-.3)
ax.set_ylim(0, 2.6)
ax.set_yticks(np.arange(0, 2.6, .5))
ax.tick_params(axis='both', length=0, pad=8, colors='#344554')
ax.yaxis.grid(True, color='#E3E8EC', linewidth=.75, zorder=0)
ax.set_ylabel('Boltz2 affinity score · log₁₀(IC50 / µM)', labelpad=13, color='#344554')
ax.set_xlabel('Amino-acid substitution at canonical residue 177', labelpad=13, color='#344554')
ax.legend(loc='upper right', bbox_to_anchor=(1, 1.12), frameon=False, ncol=2,
          handlelength=1.6, columnspacing=2)
fig.text(.078, .925, 'EPSPS G177 substitution scan', fontsize=21,
         fontweight='bold', color='#203444')
fig.text(.078, .871, 'Predicted PEP and glyphosate affinity across all 19 substitutions and wild type',
         fontsize=12, color='#465C6D')
fig.text(.078, .807, 'Lower score = stronger predicted binding',
         fontsize=11, color='#344554')
fig.text(.078, .090, 'Bars: mean of 3 prediction seeds. Error bars: sample SD. Shaded group: wild-type glycine. S3P included in both contexts.',
         fontsize=9.5, color='#465C6D')
fig.text(.078, .049, 'Exploratory Boltz 2.2.1 scores, not measured Kd or evidence of catalytic activity. Affinity-scored ligands were neutralized by the parser.',
         fontsize=9.5, color='#465C6D')
for extension in ('png', 'svg'):
    target = args.output_dir / f'g177-affinity-column-plot.{extension}'
    fig.savefig(target, dpi=220)
    print(target)
plt.close(fig)
