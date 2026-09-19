from __future__ import annotations

import csv
import io

from ..backends.interfaces import FunctionRetentionBackend
from ..schemas.models import (
    FunctionRetentionRecord, Ligand, MutationCandidate, Provenance, ScorePacket, StructureModel, TargetProtein,
)


REQUIRED_STRUCTURE_METRICS = (
    "query_tm_score", "target_tm_score", "alignment_lddt", "alignment_coverage",
    "active_site_rmsd_angstrom",
)


class FunctionRetentionAgent:
    """Evaluate negative design while keeping structure, affinity, and proxy scores separate."""

    def __init__(
        self,
        backend: FunctionRetentionBackend | None = None,
        min_tm_score: float = 0.80,
        min_alignment_lddt: float = 0.70,
        min_alignment_coverage: float = 0.80,
        max_active_site_rmsd_angstrom: float = 1.50,
        min_herbicide_kd_fold_change: float = 10.0,
        max_native_kd_fold_change: float = 3.0,
        max_fold_ddg_kcal_mol: float = 2.0,
    ) -> None:
        self.backend = backend
        self.min_tm_score = min_tm_score
        self.min_alignment_lddt = min_alignment_lddt
        self.min_alignment_coverage = min_alignment_coverage
        self.max_active_site_rmsd_angstrom = max_active_site_rmsd_angstrom
        self.min_herbicide_kd_fold_change = min_herbicide_kd_fold_change
        self.max_native_kd_fold_change = max_native_kd_fold_change
        self.max_fold_ddg_kcal_mol = max_fold_ddg_kcal_mol

    def assess(
        self,
        target: TargetProtein,
        candidates: list[MutationCandidate],
        scores: list[ScorePacket],
        herbicide: Ligand,
        native_ligands: list[Ligand],
        wild_type_structures: list[StructureModel],
    ) -> list[FunctionRetentionRecord]:
        records = []
        for candidate, score in zip(candidates, scores):
            evidence = self.backend.evaluate(
                target, candidate, herbicide, native_ligands, wild_type_structures
            ) if self.backend else {}
            structural = {name: evidence.get("structural_metrics", {}).get(name) for name in REQUIRED_STRUCTURE_METRICS}
            kd = {herbicide.name: evidence.get("ligand_kd_molar", {}).get(herbicide.name)}
            fold_change = {
                herbicide.name: evidence.get("ligand_kd_fold_change_vs_wt", {}).get(herbicide.name)
            }
            for ligand in native_ligands:
                kd[ligand.name] = evidence.get("ligand_kd_molar", {}).get(ligand.name)
                fold_change[ligand.name] = evidence.get("ligand_kd_fold_change_vs_wt", {}).get(ligand.name)
            ddg = evidence.get("fold_ddg_kcal_mol")
            missing = [f"structural_metrics.{name}" for name, value in structural.items() if value is None]
            missing += [f"ligand_kd_molar.{name}" for name, value in kd.items() if value is None]
            missing += [f"ligand_kd_fold_change_vs_wt.{name}" for name, value in fold_change.items() if value is None]
            if ddg is None:
                missing.append("fold_ddg_kcal_mol")
            decision = "INSUFFICIENT_EVIDENCE" if missing else self._decision(
                structural, fold_change, herbicide.name, [ligand.name for ligand in native_ligands], float(ddg)
            )
            records.append(FunctionRetentionRecord(
                candidate.mutation, structural, kd, fold_change, ddg,
                score.herbicide_escape_score, score.native_ligand_retention_score,
                decision, missing,
                list(evidence.get("provenance", [])) or [Provenance(
                    "computed://function-retention-gate", "evidence-completeness-gate", "evaluation",
                    "No Kd or mutant structural evidence was inferred from heuristic scores",
                )],
            ))
        return records

    def _decision(self, structural, fold_change, herbicide, native_ligands, ddg):
        structure_ok = (
            structural["query_tm_score"] >= self.min_tm_score
            and structural["target_tm_score"] >= self.min_tm_score
            and structural["alignment_lddt"] >= self.min_alignment_lddt
            and structural["alignment_coverage"] >= self.min_alignment_coverage
            and structural["active_site_rmsd_angstrom"] <= self.max_active_site_rmsd_angstrom
        )
        affinity_ok = (
            fold_change[herbicide] >= self.min_herbicide_kd_fold_change
            and all(fold_change[name] <= self.max_native_kd_fold_change for name in native_ligands)
        )
        return (
            "MEETS_COMPUTATIONAL_SCREEN"
            if structure_ok and affinity_ok and ddg <= self.max_fold_ddg_kcal_mol
            else "FAILS_COMPUTATIONAL_SCREEN"
        )


def retention_csv(records: list[FunctionRetentionRecord], herbicide: str, native_ligands: list[str]) -> str:
    output = io.StringIO()
    fields = ["mutation", "decision", "query_tm_score", "target_tm_score", "alignment_lddt", "alignment_coverage",
              "active_site_rmsd_angstrom",
              f"{herbicide}_kd_molar", f"{herbicide}_kd_fold_change_vs_wt"]
    fields += [item for ligand in native_ligands for item in (f"{ligand}_kd_molar", f"{ligand}_kd_fold_change_vs_wt")]
    fields += ["fold_ddg_kcal_mol", "herbicide_escape_score", "native_ligand_retention_score", "missing_evidence"]
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for record in records:
        row = {
            "mutation": record.mutation, "decision": record.decision,
            **record.structural_metrics,
            f"{herbicide}_kd_molar": record.ligand_kd_molar.get(herbicide),
            f"{herbicide}_kd_fold_change_vs_wt": record.ligand_kd_fold_change_vs_wt.get(herbicide),
            "fold_ddg_kcal_mol": record.fold_ddg_kcal_mol,
            "herbicide_escape_score": record.herbicide_escape_score,
            "native_ligand_retention_score": record.native_ligand_retention_score,
            "missing_evidence": ";".join(record.missing_evidence),
        }
        for ligand in native_ligands:
            row[f"{ligand}_kd_molar"] = record.ligand_kd_molar.get(ligand)
            row[f"{ligand}_kd_fold_change_vs_wt"] = record.ligand_kd_fold_change_vs_wt.get(ligand)
        writer.writerow(row)
    return output.getvalue()


def retention_markdown(records: list[FunctionRetentionRecord], herbicide: str, native_ligands: list[str]) -> str:
    headings = ["Mutation", "Decision", "qTM", "LDDT", "Active-site RMSD (Å)", f"{herbicide} Kd"]
    headings += [f"{name} Kd" for name in native_ligands]
    headings += [f"{herbicide} Kd/WT", "Missing evidence"]
    lines = ["| " + " | ".join(headings) + " |", "|" + "|".join(["---"] * len(headings)) + "|"]
    for record in records:
        values = [
            record.mutation, record.decision,
            _display(record.structural_metrics.get("query_tm_score")),
            _display(record.structural_metrics.get("alignment_lddt")),
            _display(record.structural_metrics.get("active_site_rmsd_angstrom")),
            _display_kd(record.ligand_kd_molar.get(herbicide)),
        ]
        values += [_display_kd(record.ligand_kd_molar.get(name)) for name in native_ligands]
        values += [_display(record.ligand_kd_fold_change_vs_wt.get(herbicide)), ", ".join(record.missing_evidence)]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines) + "\n"


def _display(value: float | None) -> str:
    return "N/A" if value is None else f"{value:.4g}"


def _display_kd(value: float | None) -> str:
    return "N/A" if value is None else f"{value * 1e6:.4g} µM"
