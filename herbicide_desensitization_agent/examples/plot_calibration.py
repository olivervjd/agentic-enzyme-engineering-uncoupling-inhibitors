"""Render calibration variability and experimental-reference errors with a full legend."""
import argparse
import json
import textwrap
from pathlib import Path


def render(root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = json.loads((root / "calibration_report.json").read_text())
    selections = [
        ("Arabidopsis WT repeatability", [r for r in report["comparisons"]
                                         if r["subject"] == "arabidopsis_wt" and r["kind"] == "repeatability"]),
        ("E. coli experimental-reference errors", [r for r in report["comparisons"]
                                                  if r["kind"] == "experimental_accuracy_homolog"]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 7.2))
    labels = {"herbicide": "Glyphosate + S3P", "native": "PEP + S3P", "native_s3p": "PEP + S3P\nS3P scored"}
    for ax, (title, rows) in zip(axes, selections):
        groups = sorted({(r["subject"], r["context"]) for r in rows})
        for index, group in enumerate(groups):
            subset = [r for r in rows if (r["subject"], r["context"]) == group]
            for key, offset, color, label in (("global_ca_rmsd_angstrom", -.12, "#176b87", "Global CA"),
                                               ("active_site_rmsd_angstrom", .12, "#b04a56", "Pocket backbone")):
                values = [r[key] for r in subset if r[key] is not None]
                ax.scatter([index + offset + .035 * (i - (len(values) - 1) / 2) for i in range(len(values))],
                           values, color=color, s=40, alpha=.8, label=label if index == 0 else None)
        ax.axhline(1.0, linestyle="--", linewidth=1, color="#666666", label="Original 1 A limit")
        ax.set_xticks(range(len(groups)), [labels[c] if s == "arabidopsis_wt" else
                                          ("WT / 1G6S\nGlyphosate + S3P" if s == "ecoli_wt" else "G96A / 1MI4\nS3P only")
                                          for s, c in groups])
        ax.set_title(title, fontsize=12)
        ax.set_ylabel("RMSD (angstroms)")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=.2)
        ax.legend(fontsize=8)
    caption = ("Figure legend: each dot is one globally fitted structural comparison. Left: all three pairs of three "
               "seeded WT predictions in each ligand context. Right: each of three E. coli control predictions versus "
               "its matched-context experimental structure. Blue shows CA RMSD; red shows active-site N/CA/C/O RMSD "
               "using the same global fit. Lower values mean closer agreement. Dashed lines retain the original "
               "screening limit, not a calibrated error bound. All predictions use validated MSAs and no templates. "
               "These historical controls are not verified as training-set holdouts. Homolog controls do not establish "
               "Arabidopsis accuracy; no panel demonstrates retained enzyme "
               "activity or herbicide resistance. Ligand protonation and independent S3P-aware docking remain unresolved.")
    fig.suptitle("EPSPS calibration: repeatability is not accuracy", fontsize=16, y=.97)
    fig.subplots_adjust(left=.07, right=.98, top=.86, bottom=.34, wspace=.26)
    fig.text(.07, .04, textwrap.fill(caption, 145), fontsize=9, va="bottom", linespacing=1.4)
    fig.savefig(root / "calibration_structure_checks.png", dpi=180)
    fig.savefig(root / "calibration_structure_checks.pdf")
    plt.close(fig)
    (root / "calibration_figure_legend.md").write_text(caption + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    render(parser.parse_args().root)
