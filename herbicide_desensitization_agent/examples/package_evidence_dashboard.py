"""Package the built workbench as an offline HTML report and annotated figures."""
from __future__ import annotations

import argparse
import html
import json
import re
import shutil
from pathlib import Path

from ..app.backends.calibration import sha, write


def figure_scope_rows(data):
    """Keep pose-only records and independent docking out of affinity/structure plots."""
    binding = [r for r in data['tables']['binding'] if r['variant']=='arabidopsis_wt'
               and r.get('affinity_metric')=='predicted pIC50' and r.get('replicate_results')
               and all('value' in value for value in r['replicate_results'])]
    contacts = [r for r in data['tables']['residue_interactions'] if r['variant']=='arabidopsis_wt'
                and r.get('model_mode')!='docking' and r.get('method')!='AutoDock Vina'
                and isinstance(r.get('residue'),int)]
    same_target = [r for r in data['structural_comparisons'] if r.get('kind') in
                   {'experimental_accuracy_same_target','experimental_accuracy_same_target_apo'}]
    return binding, contacts, same_target


def figures(root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data = json.loads((root / "data/evidence_system.json").read_text())
    folder = root / "figures"
    folder.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none", "savefig.dpi": 300})
    def save(fig, name, legend):
        fig.savefig(folder / f"{name}.png", bbox_inches="tight", facecolor="white")
        fig.savefig(folder / f"{name}.svg", bbox_inches="tight", facecolor="white")
        (folder / f"{name}.legend.txt").write_text(legend + "\n")
        plt.close(fig)
    binding, contacts, same_target = figure_scope_rows(data)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, row in enumerate(binding):
        values = [r["value"] for r in row["replicate_results"]]
        ax.scatter([i + (j - (len(values)-1)/2) * .055 for j in range(len(values))], values, color=["#bc4841", "#3275aa", "#568e81"][i % 3], s=65)
        ax.hlines(row["estimate"], i - .18, i + .18, color="#1d3544", linewidth=2)
    ax.set_xticks(range(len(binding)), [r["ligand"].upper() for r in binding])
    ax.set_ylabel("Predicted pIC50 (dimensionless)")
    ax.set_title("WT EPSPS · separate ligand predictions", loc="left", fontweight="bold", pad=18)
    counts=', '.join(f"{r['ligand'].upper()}: n={len(r['replicate_results'])}" for r in binding)
    fig.text(.12, -.02, f"Recorded affinity seeds ({counts}); dots are predictions, bars are means.\nNo calibrated prediction interval. These values do not establish binding or native function.", fontsize=9)
    save(fig, "ligand-affinity", data["legends"]["binding"])
    positions = sorted({r["residue"] for r in contacts})
    values = np.full((len(positions), 3), np.nan)
    for r in contacts:
        if r["contact_frequency"] is not None:
            values[positions.index(r["residue"]), ["glyphosate", "pep", "s3p"].index(r["ligand"])] = r["contact_frequency"]
    fig, ax = plt.subplots(figsize=(6.5, max(6, len(positions) * .22)))
    im = ax.imshow(values, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(3), ["Glyphosate", "PEP + S3P", "S3P + PEP"])
    ax.set_yticks(range(len(positions)), positions)
    ax.set_ylabel("Canonical Arabidopsis residue")
    ax.set_title("WT Boltz structure-mode contact frequency", loc="left", fontweight="bold", pad=18)
    denominators=sorted({r['replicates'] for r in contacts if isinstance(r.get('replicates'),int)})
    fig.colorbar(im, ax=ax, label="Observed seed fraction (n="+','.join(map(str,denominators))+")", shrink=.55)
    fig.text(.06, .01, "Heavy-atom proximity ≤5 Å; no binding claim. Independent Vina contacts remain in separate tables.", fontsize=8)
    save(fig, "residue-contact-matrix", "Boltz structure-mode records only; Vina observations are not pooled. " + data["contact_definitions"])
    comparisons = data["structural_comparisons"]
    groups = [("WT seed comparisons", [r for r in comparisons if r.get("subject") == "arabidopsis_wt"
                and r.get("kind") in {"repeatability", "wt_repeatability"} and "query_tm_score" in r]),
              ("Same-target apo crystal", same_target),
              ("Homolog vs crystal", [r for r in comparisons if r.get("kind") == "experimental_accuracy_homolog"])]
    groups = [(name, rows) for name, rows in groups if rows]
    metrics = [("query_tm_score", "Directional TM-score", False), ("alignment_lddt", "Cα lDDT", False),
               ("global_ca_rmsd_angstrom", "Global Cα RMSD, Å", True), ("active_site_rmsd_angstrom", "Active-site RMSD, Å", True)]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, (metric, label, lower) in zip(axes.flat, metrics):
        for i, (name, rows) in enumerate(groups):
            vals = [r[metric] for r in rows if r.get(metric) is not None]
            if vals:
                ax.scatter([i + (j % 5 - 2) * .035 for j in range(len(vals))], vals, s=20, color="#287d86", alpha=.75)
        limit = data["calibration"]["thresholds"].get(metric)
        if limit is not None:
            ax.axhline(limit, color="#bd4e48", linestyle="--", linewidth=1, label="Historical cutoff (not calibrated)")
        ax.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=8)
        ax.set_title(label, loc="left", fontsize=11)
        if ax is axes.flat[0] and ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7, loc="best")
    fig.suptitle("Replicate agreement and experimental accuracy retain separate scopes", x=.08, ha="left", fontweight="bold")
    fig.tight_layout(rect=[0, .04, 1, .94])
    fig.text(.08, .005, "Paired comparisons are correlated. Structural agreement alone does not establish native function.", fontsize=9)
    save(fig, "structure-calibration", data["legends"]["structure"])


def package(root, dist, *, render_figures=True):
    root, dist = Path(root), Path(dist)
    index = (dist / "index.html").read_text()
    # Vite produces one self-contained ES module, no external CDN libraries.
    def script(match):
        path = dist / match.group(1).lstrip("./")
        code = path.read_text().replace("</script", "<\\/script")
        if re.search(r'from\s*["\']\./', code):
            raise ValueError("Static report requires a single JS bundle; rebuild with inlineDynamicImports")
        return '<script type="module">' + code + '</script>'
    index = re.sub(r'<script[^>]+src="/?([^" ]+)"[^>]*></script>', script, index)
    def style(match):
        return "<style>" + (dist / match.group(1).lstrip("./")).read_text() + "</style>"
    index = re.sub(r'<link[^>]+href="/?([^" ]+\.css)"[^>]*>', style, index)
    index = re.sub(r'<link[^>]+rel="modulepreload"[^>]*>', '', index)
    payload = {name: json.loads((root / "data" / filename).read_text())
               for name, filename in (("results", "results.json"), ("structures", "structures.json"), ("evidence", "evidence_system.json"))}
    review_path = root / "data/evidence_model_review.json"
    if review_path.exists():
        from .review_evidence_bundle import fingerprint
        review = json.loads(review_path.read_text())
        if (review.get("source_bundle_sha256") != sha(root / "data/evidence_system.json")
                or review.get('source_content_sha256') != fingerprint(payload['evidence'])
                or review.get('source_snapshot_current') is not True
                or review.get('status') != 'COMPLETED'):
            raise ValueError("Model commentary source changed; rerun review or omit stale commentary")
        payload["evidence"]["model_review"] = review
    payload["evidence"]["static_report"] = True
    encoded = json.dumps(payload, allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    index = index.replace("<head>", '<head><script>window.__EVIDENCE_DATA__=' + encoded + ";</script>", 1)
    (root / "report.html").write_text(index)
    if render_figures:
        figures(root)
    shutil.copytree(dist, root / "dashboard", dirs_exist_ok=True)
    shutil.copytree(root / "data", root / "dashboard/data", dirs_exist_ok=True)
    readme = """# Evidence workbench

Open `report.html` for the self-contained offline interactive report. It embeds
the JavaScript, styles, molecular geometry and all displayed result tables.
The dashboard directory can also be served with `python3 -m http.server 8765`.

The current scientific decision is **INSUFFICIENT_EVIDENCE**. No new mutations
are nominated; no computational result is labelled herbicide resistant.

CSV/JSON tables and legends are in data/. Verified input structures, MSA and
reference evidence are in raw/. Publication figures are PNG and editable SVG.
workflow-manifest.json records source hashes, methods, versions and artifact hashes.
Read limitations.md and experimental-validation-plan.md before interpretation.
"""
    (root / "README.md").write_text(readme)
    manifest_path = root / "workflow-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["artifacts"] = [{"path": str(p.relative_to(root)), "sha256": sha(p)} for p in sorted(root.rglob("*"))
                             if p.is_file() and p != manifest_path]
    write(manifest_path, manifest)
    print("Packaged offline report, dashboard and three PNG/SVG figure pairs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dist", type=Path, required=True)
    args = parser.parse_args()
    package(args.output, args.dist)


if __name__ == "__main__":
    main()
