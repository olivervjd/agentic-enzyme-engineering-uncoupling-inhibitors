"""Render campaign structural comparisons with WT variability and a complete figure legend."""
import argparse
import json
import textwrap
from pathlib import Path

from .evaluate_epsps_campaign import COMPARISON_LEGEND


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = json.loads((args.results / "structure_comparisons.json").read_text())
    campaign = json.loads((args.results / "campaign_manifest.json").read_text())
    thresholds = json.loads((args.results / "function_retention_thresholds.json").read_text())
    labels = ["WT-control"] + list(dict.fromkeys(row["mutation"] for row in rows if row["mutation"] not in {"WT", "WT-control"}))
    metrics = [
        ("query_tm_score", "A  TM-align score", thresholds["min_tm_score"], "higher", (0.88, 1.01)),
        ("alignment_lddt", "B  CA lDDT", thresholds["min_alignment_lddt"], "higher", (0.70, 1.01)),
        ("global_ca_rmsd_angstrom", "C  Global CA RMSD (A)", thresholds["max_global_ca_rmsd_angstrom"], "lower", None),
        ("active_site_rmsd_angstrom", "D  Active-site backbone RMSD (A)", thresholds["max_active_site_rmsd_angstrom"], "lower", None),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.85, bottom=0.31, hspace=0.40, wspace=0.20)
    conditions = campaign.get("conditions", ["native", "herbicide"])
    colors = {"native": "#007f79", "herbicide": "#bf4162", "native_s3p": "#8c6b18"}
    names_by_condition = {"native": "PEP + S3P", "herbicide": "Glyphosate + S3P", "native_s3p": "S3P affinity context"}
    for ax, (key, title, threshold, direction, limits) in zip(axes.flat, metrics):
        for index, label in enumerate(labels):
            for ci, condition in enumerate(conditions):
                offset = (ci - (len(conditions) - 1) / 2) * .22
                values = [row[key] for row in rows if row["mutation"] == label and row["condition"] == condition]
                spread = [offset + (i - (len(values) - 1) / 2) * 0.045 for i in range(len(values))]
                ax.scatter([index + value for value in spread], values, s=25, color=colors[condition],
                           alpha=0.85, edgecolor="white", linewidth=0.4,
                           label=names_by_condition[condition] if index == 0 else None)
        ax.axhline(threshold, color="#555555", linewidth=1, linestyle="--", label=f"Screen threshold ({direction} is better)")
        ax.set_title(title, loc="left", fontsize=12)
        ax.set_xticks(range(len(labels)), ["WT variability"] + labels[1:], rotation=25, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        if limits:
            observed = [row[key] for row in rows if row["mutation"] != "WT"]
            ax.set_ylim(min(limits[0], min(observed) - 0.02), limits[1])
        else:
            ax.set_ylim(bottom=0)
    fig.suptitle("EPSPS mutant structure screen", x=0.07, ha="left", fontsize=20, fontweight="bold")
    replicates = len([r for r in campaign["records"] if r["mutation"] == "WT" and r["condition"] == "native"])
    fig.text(0.07, 0.93, f"{len(campaign['records'])} Boltz-2 models | {len(labels) - 1} substitutions | "
             f"{len(conditions)} contexts | {replicates} runs per context", fontsize=11)
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles[:len(conditions)], names[:len(conditions)], loc="upper right", bbox_to_anchor=(0.98, 0.905), frameon=False, ncol=len(conditions))
    legend = "Figure 1. " + COMPARISON_LEGEND + " Dashed lines show the predeclared computational thresholds."
    fig.text(0.07, 0.235, textwrap.fill(legend, 158), va="top", fontsize=9, linespacing=1.45)
    fig.savefig(args.results / "structural_comparison.png", dpi=180)
    fig.savefig(args.results / "structural_comparison.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
