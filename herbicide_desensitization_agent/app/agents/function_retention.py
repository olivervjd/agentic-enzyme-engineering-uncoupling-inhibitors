from __future__ import annotations

import csv
import io
import math

from ..backends.interfaces import FunctionRetentionBackend
from ..schemas.models import FunctionRetentionRecord, Provenance


REQUIRED_STRUCTURE_METRICS = (
    "query_tm_score", "target_tm_score", "alignment_lddt", "alignment_coverage",
    "active_site_rmsd_angstrom", "global_ca_rmsd_angstrom",
)

RETENTION_LEGEND = (
    "Table legend: WT is the unmutated reference, not a resistant candidate. Mutation labels use "
    "the supplied target sequence numbering. Kd is a dissociation constant in mol/L (M); larger "
    "values mean weaker binding. Bounds are supplied prediction intervals, with their interpretation "
    "recorded per row; they are not automatically confidence intervals. Kd/WT is computed from "
    "matched estimates, never from docking confidence, heuristic scores, or IC50. Ratio bounds are "
    "[mutant lower / WT upper, mutant upper / WT lower]. Native binding passes only when the entire "
    "ratio interval lies inside the configured two-sided equivalence window for EVERY native ligand. "
    "Herbicide escape requires its entire interval above the configured weakening threshold. "
    "qTM and tTM are TM-align scores normalized by mutant and reference lengths (0-1; higher is better). "
    "CA lDDT measures local distance agreement (0-1); coverage is the minimum fraction of each modeled "
    "chain matched. Global CA RMSD and active-site backbone RMSD are in angstroms after global "
    "superposition (lower is better). Metrics concern modeled residues, not unmodeled segments. "
    "WT self-comparison is an identity control, not experimental validation. Fold ddG is mutant minus "
    "WT folding free energy in kcal/mol. Missing data are N/A (blank in CSV), never zero. A computational "
    "pass does not establish absence of herbicide binding, identical structure, preserved catalytic "
    "turnover, signaling, or herbicide resistance. Synthetic fixtures are software tests only."
)


def _number(value, name, *, positive=False, unit_interval=False):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    if unit_interval and not 0 <= value <= 1:
        raise ValueError(f"{name} must be in [0, 1]")
    return float(value)


class FunctionRetentionAgent:
    """Require matched WT evidence and uncertainty-aware equivalence, independently of proxy scores."""

    def __init__(
        self, backend: FunctionRetentionBackend | None = None,
        min_tm_score=0.95, min_alignment_lddt=0.90, min_alignment_coverage=0.95,
        max_active_site_rmsd_angstrom=1.0, min_herbicide_kd_fold_change=10.0,
        max_native_kd_fold_change=3.0, max_fold_ddg_kcal_mol=2.0,
        max_global_ca_rmsd_angstrom=1.0,
    ) -> None:
        self.backend = backend
        self.thresholds = {
            "min_tm_score": min_tm_score, "min_alignment_lddt": min_alignment_lddt,
            "min_alignment_coverage": min_alignment_coverage,
            "max_active_site_rmsd_angstrom": max_active_site_rmsd_angstrom,
            "max_global_ca_rmsd_angstrom": max_global_ca_rmsd_angstrom,
            "min_herbicide_kd_fold_change": min_herbicide_kd_fold_change,
            "max_native_kd_fold_change": max_native_kd_fold_change,
            "max_fold_ddg_kcal_mol": max_fold_ddg_kcal_mol,
        }
        for name, value in self.thresholds.items():
            _number(value, name, positive=True,
                    unit_interval=name.startswith("min_") and "fold_change" not in name)
        if max_native_kd_fold_change < 1 or min_herbicide_kd_fold_change <= 1:
            raise ValueError("Affinity thresholds must define equivalence and weaker herbicide binding")

    def _read(self, evidence, ligands):
        structural = {}
        for name in REQUIRED_STRUCTURE_METRICS:
            value = _number(evidence.get("structural_metrics", {}).get(name), name,
                            unit_interval="rmsd" not in name)
            if value is not None and value < 0:
                raise ValueError(f"{name} must be nonnegative")
            structural[name] = value
        kd, intervals = {}, {}
        for name in ligands:
            kd[name] = _number(evidence.get("ligand_kd_molar", {}).get(name), f"Kd/{name}", positive=True)
            bounds = evidence.get("ligand_kd_intervals_molar", {}).get(name)
            if bounds is not None:
                if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
                    raise ValueError(f"Expected two interval bounds for {name}")
                bounds = [_number(x, f"interval/{name}", positive=True) for x in bounds]
                if any(x is None for x in bounds) or bounds[0] > bounds[1]:
                    raise ValueError(f"Invalid interval for {name}")
                if kd[name] is None or not bounds[0] <= kd[name] <= bounds[1]:
                    raise ValueError(f"Kd must lie in its interval for {name}")
            intervals[name] = bounds
        ddg = _number(evidence.get("fold_ddg_kcal_mol"), "fold_ddg_kcal_mol")
        return structural, kd, intervals, ddg

    def assess(self, target, candidates, scores, herbicide, native_ligands, wild_type_structures):
        if len(candidates) != len(scores):
            raise ValueError("Each candidate requires exactly one score packet")
        if len({item.mutation for item in candidates}) != len(candidates):
            raise ValueError("Duplicate mutation identity")
        names = [herbicide.name] + [item.name for item in native_ligands]
        if not native_ligands or len(names) != len(set(names)):
            raise ValueError("Distinct herbicide and native ligand identities are required")
        wt = self.backend.reference(target, herbicide, native_ligands, wild_type_structures) if self.backend else {}
        wt_structure, wt_kd, wt_intervals, _ = self._read(wt, names)
        fallback = [Provenance(
            "computed://function-retention-gate", "evidence-completeness-gate", "evaluation",
            "Unavailable measurements are not inferred from heuristic scores",
        )]
        wt_missing = [f"WT.ligand_kd_molar.{name}" for name in names if wt_kd[name] is None]
        wt_missing += [f"WT.ligand_kd_intervals_molar.{name}" for name in names if wt_intervals[name] is None]
        wt_missing += [f"WT.structural_metrics.{key}" for key, value in wt_structure.items() if value is None]
        records = [FunctionRetentionRecord(
            "WT", wt_structure, wt_kd, {name: 1.0 if wt_kd[name] else None for name in names},
            0.0, 0.0, 0.0, "WILD_TYPE_REFERENCE", wt_missing, wt.get("provenance") or fallback,
            wt_intervals, {name: [1.0, 1.0] if wt_intervals[name] else None for name in names},
            affinity_protocol=wt.get("affinity_protocol", {}),
            interval_description=wt.get("interval_description", ""),
            structural_method=wt.get("structural_method", "unavailable"),
        )]
        for candidate, score in zip(candidates, scores):
            evidence = self.backend.evaluate(target, candidate, herbicide, native_ligands, wild_type_structures) if self.backend else {}
            structural, kd, intervals, ddg = self._read(evidence, names)
            ratios, ratio_intervals = {}, {}
            missing = [f"structural_metrics.{name}" for name, value in structural.items() if value is None]
            if not evidence.get("structural_method"):
                missing.append("structural_method")
            for label, item in (("WT", wt), (candidate.mutation, evidence)):
                for field in ("provenance", "interval_description"):
                    if not item.get(field):
                        missing.append(f"{label}.{field}")
            for name in names:
                ratios[name] = kd[name] / wt_kd[name] if kd[name] is not None and wt_kd[name] is not None else None
                ratio_intervals[name] = (
                    [intervals[name][0] / wt_intervals[name][1], intervals[name][1] / wt_intervals[name][0]]
                    if intervals[name] and wt_intervals[name] else None
                )
                for label, values in (("WT.ligand_kd_molar", wt_kd), ("ligand_kd_molar", kd),
                                      ("WT.ligand_kd_intervals_molar", wt_intervals), ("ligand_kd_intervals_molar", intervals)):
                    if values[name] is None:
                        missing.append(f"{label}.{name}")
                protocol = evidence.get("affinity_protocol", {}).get(name)
                if not protocol or protocol != wt.get("affinity_protocol", {}).get(name):
                    missing.append(f"matched_affinity_protocol.{name}")
                supplied = evidence.get("ligand_kd_fold_change_vs_wt", {}).get(name)
                if supplied is not None:
                    _number(supplied, f"supplied ratio/{name}", positive=True)
                    if ratios[name] is not None and not math.isclose(supplied, ratios[name], rel_tol=1e-6):
                        raise ValueError(f"Supplied fold change disagrees with WT Kd for {candidate.mutation}/{name}")
            if ddg is None:
                missing.append("fold_ddg_kcal_mol")
            failed = self._failed_checks(structural, ratio_intervals, herbicide.name, names[1:], ddg)
            decision = "FAILS_COMPUTATIONAL_SCREEN" if failed else (
                "INSUFFICIENT_EVIDENCE" if missing else "MEETS_COMPUTATIONAL_SCREEN"
            )
            records.append(FunctionRetentionRecord(
                candidate.mutation, structural, kd, ratios, ddg,
                score.herbicide_escape_score, score.native_ligand_retention_score, decision, missing,
                evidence.get("provenance") or fallback, intervals, ratio_intervals, failed,
                evidence.get("affinity_protocol", {}), evidence.get("interval_description", ""),
                evidence.get("structural_method", "unavailable"),
            ))
        return records

    def _failed_checks(self, structural, intervals, herbicide, native_ligands, ddg):
        t = self.thresholds
        minimums = {"query_tm_score": t["min_tm_score"], "target_tm_score": t["min_tm_score"],
                    "alignment_lddt": t["min_alignment_lddt"], "alignment_coverage": t["min_alignment_coverage"]}
        maximums = {"active_site_rmsd_angstrom": t["max_active_site_rmsd_angstrom"],
                    "global_ca_rmsd_angstrom": t["max_global_ca_rmsd_angstrom"]}
        failed = [name for name, bound in minimums.items() if structural[name] is not None and structural[name] < bound]
        failed += [name for name, bound in maximums.items() if structural[name] is not None and structural[name] > bound]
        if intervals[herbicide] and intervals[herbicide][0] < t["min_herbicide_kd_fold_change"]:
            failed.append(f"herbicide_weakening.{herbicide}")
        for name in native_ligands:
            if intervals[name] and not (1 / t["max_native_kd_fold_change"] <= intervals[name][0]
                                       <= intervals[name][1] <= t["max_native_kd_fold_change"]):
                failed.append(f"native_equivalence.{name}")
        if ddg is not None and ddg > t["max_fold_ddg_kcal_mol"]:
            failed.append("fold_ddg_kcal_mol")
        return failed


def retention_csv(records, herbicide, native_ligands):
    output = io.StringIO()
    fields = ["mutation", "decision", "structural_method", *REQUIRED_STRUCTURE_METRICS]
    names = [herbicide] + native_ligands
    suffixes = ("kd_molar", "kd_lower_molar", "kd_upper_molar", "kd_fold_change_vs_wt", "ratio_lower", "ratio_upper", "protocol")
    for ligand in names:
        fields += [f"{ligand}_{suffix}" for suffix in suffixes]
    fields += ["fold_ddg_kcal_mol", "interval_description", "failed_checks", "missing_evidence", "provenance"]
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for record in records:
        row = {"mutation": record.mutation, "decision": record.decision, "structural_method": record.structural_method,
               **record.structural_metrics, "fold_ddg_kcal_mol": record.fold_ddg_kcal_mol,
               "interval_description": record.interval_description, "failed_checks": ";".join(record.failed_checks),
               "missing_evidence": ";".join(record.missing_evidence),
               "provenance": ";".join(f"{p.source} ({p.method}; {p.evidence_type})" for p in record.provenance)}
        for name in names:
            bounds = record.ligand_kd_intervals_molar.get(name) or [None, None]
            ratio = record.ligand_kd_ratio_intervals.get(name) or [None, None]
            values = (record.ligand_kd_molar.get(name), *bounds, record.ligand_kd_fold_change_vs_wt.get(name),
                      *ratio, record.affinity_protocol.get(name, ""))
            row.update({f"{name}_{suffix}": value for suffix, value in zip(suffixes, values)})
        writer.writerow(row)
    return output.getvalue()


def retention_markdown(records, herbicide, native_ligands):
    headings = ["Mutation", "Decision", "qTM / tTM", "CA lDDT", "Coverage", "Global / site RMSD (A)"]
    names = [herbicide] + native_ligands
    headings += [f"{name} Kd (M) [bounds]; ratio [bounds]" for name in names]
    headings += ["Failed checks", "Missing evidence"]
    lines = ["| " + " | ".join(headings) + " |", "|" + "|".join(["---"] * len(headings)) + "|"]
    for record in records:
        s = record.structural_metrics
        values = [record.mutation, record.decision,
                  f"{_display(s.get('query_tm_score'))} / {_display(s.get('target_tm_score'))}",
                  _display(s.get("alignment_lddt")), _display(s.get("alignment_coverage")),
                  f"{_display(s.get('global_ca_rmsd_angstrom'))} / {_display(s.get('active_site_rmsd_angstrom'))}"]
        for name in names:
            bounds = record.ligand_kd_intervals_molar.get(name) or [None, None]
            ratio = record.ligand_kd_ratio_intervals.get(name) or [None, None]
            values.append(f"{_display(record.ligand_kd_molar.get(name))} [{_display(bounds[0])}, {_display(bounds[1])}]; "
                          f"{_display(record.ligand_kd_fold_change_vs_wt.get(name))} [{_display(ratio[0])}, {_display(ratio[1])}]")
        values += [", ".join(record.failed_checks) or "none", ", ".join(record.missing_evidence) or "none"]
        lines.append("| " + " | ".join(str(value).replace("|", "\\|").replace("\n", " ") for value in values) + " |")
    lines += ["", RETENTION_LEGEND, "", "Methods and provenance:", ""]
    for record in records:
        lines.append(f"- {record.mutation}: structure={record.structural_method}; "
                     f"affinity protocols={record.affinity_protocol or 'unavailable'}; "
                     f"intervals={record.interval_description or 'unavailable'}; "
                     f"sources={'; '.join(p.source + ' (' + p.evidence_type + ')' for p in record.provenance)}")
    return "\n".join(lines) + "\n"


def _display(value):
    return "N/A" if value is None else f"{value:.4g}"
