"""Plot separately seeded affinity predictions without implying statistical equivalence."""
import argparse
import json
from pathlib import Path
import textwrap

FIGURE_LEGEND = (
    "Figure 2. Boltz-2 predicted pIC50 for WT and six contact-selected EPSPS substitutions. "
    "PEP is phosphoenolpyruvate and S3P is shikimate 3-phosphate. Glyphosate and PEP are each predicted "
    "with S3P present; S3P affinity is predicted with PEP present. Each point is one model prediction "
    "at seed 43 or 44, not an experimental measurement or biological replicate. "
    "Predicted pIC50 = 6 - affinity_pred_value is dimensionless; higher values mean stronger predicted "
    "affinity. These are not Kd or Km measurements. No confidence intervals or equivalence bounds are "
    "estimated from two runs. WT is the reference, not a resistant candidate. Query-only MSA and "
    "unenumerated protonation states limit interpretation. Native-function retention and herbicide "
    "resistance cannot be established by these plots."
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    root = parser.parse_args().root
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = json.loads((root / "affinity_predictions.json").read_text())["records"]
    mutations = json.loads((root / "inputs/input_manifest.json").read_text())["mutations"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 7.5), sharey=True)
    fig.subplots_adjust(left=.065, right=.98, bottom=.36, top=.83, wspace=.12)
    for axis, ligand, title in zip(axes, ("glyphosate", "pep", "s3p"),
                                  ("A  Glyphosate + S3P", "B  PEP + S3P", "C  S3P + PEP")):
        for replicate, color, marker, offset in ((1, "#007f79", "o", -.12), (2, "#bf4162", "^", .12)):
            selected = [r for r in rows if r["ligand"] == ligand and r["replicate"] == replicate]
            values = {r["mutation"]: r["predicted_pIC50"] for r in selected}
            if len(selected) != len(mutations) or set(values) != set(mutations):
                raise ValueError("Missing or duplicate plotted prediction")
            axis.scatter([i + offset for i in range(len(mutations))], [values[m] for m in mutations],
                         color=color, marker=marker, s=45, label=f"Seed {42 + replicate}", zorder=3)
        axis.axvspan(-.4, .4, color="#ececec", zorder=0)
        axis.set_title(title, loc="left", fontsize=12)
        axis.set_xticks(range(len(mutations)), mutations, rotation=35, ha="right", fontsize=10)
        axis.grid(axis="y", alpha=.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Predicted pIC50 (dimensionless)", fontsize=11)
    fig.suptitle("EPSPS predicted affinity", x=.065, ha="left", fontsize=20, fontweight="bold")
    fig.text(.065, .90, "WT controls and contact-selected mutants | two seeded predictions per ligand", fontsize=11)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(.98, .94), ncol=2, frameon=False)
    fig.text(.065, .23, textwrap.fill(FIGURE_LEGEND, 175), va="top", fontsize=9, linespacing=1.5)
    fig.savefig(root / "affinity_predictions.png", dpi=180)
    fig.savefig(root / "affinity_predictions.pdf")
    plt.close(fig)
    (root / "affinity_figure_legend.md").write_text(FIGURE_LEGEND + "\n")


if __name__ == "__main__":
    main()
