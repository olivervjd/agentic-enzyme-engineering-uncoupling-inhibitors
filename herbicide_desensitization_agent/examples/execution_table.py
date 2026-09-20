"""Render execution facts without promoting imported evidence into new work.

The table is a projection of an actual execution journal. In its absence,
recorded scientific artifacts remain imported evidence, not completed agents.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import hashlib
import json


AGENT_SPECS = (
    ("registry", "Registry", "Define sequence, ligands, native function and required context"),
    ("literature", "Literature evidence", "Curate experimental structures, controls and source measurements"),
    ("chemical_state", "Chemical state", "Validate identity, stereochemistry, charge and protonation"),
    ("structure", "Structure", "Generate and assess replicated structures"),
    ("docking", "Docking", "Generate replicated ligand poses and pose clusters"),
    ("affinity", "Affinity", "Estimate ligand-specific quantities with uncertainty"),
    ("contact_mapping", "Contact mapping", "Calculate residue and atom interaction evidence"),
    ("calibration", "Calibration", "Assess WT variability and known controls"),
    ("mutation_design", "Mutation design", "Nominate supported conservative substitutions"),
    ("stability", "Stability", "Assess folding stability and assembly integrity"),
    ("function_retention", "Function retention", "Evaluate native-function and herbicide gates"),
    ("evidence_review", "Evidence review", "Check claims against provenance and scientific evidence"),
    ("dashboard", "Dashboard", "Present evidence without changing scientific conclusions"),
    ("learning_loop", "Learning loop", "Incorporate measured assays under validation gates"),
)


def _execution_status(records):
    """A stage is succeeded only if every recorded activity actually completed."""
    states = [r["execution_status"] for r in records]
    if "FAILED" in states:
        return "failed"
    if "RUNNING" in states:
        return "running"
    if any(s.startswith("SKIPPED_") for s in states):
        return "skipped"
    if states and all(s == "COMPLETED" for s in states):
        return "succeeded"
    return "not_started"


def _evidence_status(records):
    states = {r["evidence_status"] for r in records}
    for state in ("FAILED_VALIDATION", "INSUFFICIENT_EVIDENCE", "NOT_ASSESSED", "AVAILABLE", "VALIDATED"):
        if state in states:
            return state
    return "NOT_ASSESSED"


def execution_table(journal=None, *, journal_path=None, imported_sources=()):
    """Return fourteen explicit execution/evidence records plus run provenance.

    `journal` must conform to the shared journal validator. Unknown stage IDs are
    not heuristically assigned to an agent; they remain in the execution audit.
    `imported_sources` are existing artifacts, never current execution evidence.
    """
    from ..app.orchestrator.execution import validate_journal

    if journal is not None:
        validate_journal(journal)
    events = journal.get("stages", []) if journal else []
    provenance = None
    if journal_path:
        path = Path(journal_path)
        provenance = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    imported = list(imported_sources)
    rows = []
    for stage, name, purpose in AGENT_SPECS:
        records = [r for r in events if r["stage"] == stage]
        modes = sorted({r["execution_mode"] for r in records})
        measured = [r["runtime_seconds"] for r in records if r["execution_mode"] == "executed"]
        # Running jobs and missing timings have no final aggregate runtime.
        runtime = sum(measured) if measured and all(v is not None for v in measured) else None
        status = _execution_status(records)
        warnings = [r["reason"] for r in records if r.get("reason")]
        if not records:
            warnings.append("No current execution event recorded; existing scientific outputs do not prove this agent ran.")
        elif "imported" in modes:
            warnings.append("Imported artifacts are prior evidence, not a new model or scientific-tool invocation.")
        rows.append({
            "agent_id": stage, "agent": name, "purpose": purpose,
            "execution_status": status, "evidence_status": _evidence_status(records),
            "execution_mode": modes[0] if len(modes) == 1 else "mixed" if modes else "not_started",
            "journal_execution_status": [r["execution_status"] for r in records],
            "inputs": [a for r in records for a in r.get("input_artifacts", [])],
            "outputs": [a for r in records for a in r.get("output_artifacts", [])],
            "tools_models": [r.get("tools_models") or r["tools"] for r in records if r.get("tools_models") or r.get("tools")],
            "runtime": runtime, "runtime_unit": "seconds", "warnings": warnings,
            "provenance": {"journal": provenance, "run_id": journal.get("run_id") if journal else None,
                           "event_ids": [r["id"] for r in records]},
            "execution_records": records,
        })
    known = {stage for stage, _, _ in AGENT_SPECS}
    summary = {"run_id": journal.get("run_id") if journal else None,
               "journal_status": journal.get("status") if journal else "NOT_STARTED",
               "journal": provenance,
               "counts": dict(Counter(r["execution_status"] for r in rows)),
               "imported_evidence_sources": imported,
               "unmapped_events": [r for r in events if r["stage"] not in known],
               "claim": "Execution success and scientific evidence validity are separate; imports are not new computations."}
    return rows, summary
