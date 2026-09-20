"""Measured, durable execution events, separate from scientific evidence verdicts.

An import is not a model invocation and a skipped stage is not completed work.
JSON payload digests describe the actual supplied/returned values; file digests
describe bytes present during execution. No timestamps are inferred from files.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
import uuid
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def json_value(value):
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in sorted(value.items(), key=lambda item: str(item[0]))}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted((json_value(v) for v in value), key=lambda v: json.dumps(v, sort_keys=True))
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Execution inputs/outputs must be explicit serializable data, not {type(value).__name__}")


def payload_hash(value):
    return hashlib.sha256(json.dumps(json_value(value), allow_nan=False, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def hash_artifacts(paths):
    output = []
    for raw in paths:
        path = Path(raw).resolve()
        if not path.is_file():
            raise FileNotFoundError("Declared execution artifact does not exist")
        sha = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                sha.update(chunk)
        output.append({"path": str(path), "sha256": sha.hexdigest(), "size_bytes": path.stat().st_size})
    return output


def validate_journal(value):
    """Reject malformed execution claims before export; does not authenticate authors."""
    if not isinstance(value, dict) or value.get("schema_version") != "1.0" or not isinstance(value.get("stages"), list):
        raise ValueError("Invalid execution journal schema")
    statuses = {"RUNNING", "COMPLETED", "FAILED", "IMPORTED"} | ExecutionJournal.SKIP_STATES
    identities = set()
    def timestamp(text):
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            raise ValueError("Journal timestamps must have timezone information")
        return parsed
    timestamp(value["started_at"])
    for row in value["stages"]:
        if not isinstance(row, dict) or not row.get("id") or row["id"] in identities or not row.get("stage"):
            raise ValueError("Unique stage event identities required")
        identities.add(row["id"])
        status = row.get("execution_status")
        if status not in statuses or row.get("evidence_status") not in ExecutionJournal.EVIDENCE_STATES:
            raise ValueError("Unrecognized execution/evidence status")
        for field in ("inputs_sha256", "outputs_sha256"):
            digest = row.get(field)
            if digest is not None and (not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)):
                raise ValueError("Invalid canonical payload digest")
        if not row.get("inputs_sha256"):
            raise ValueError("Every event requires an input payload digest")
        if status in ("COMPLETED", "FAILED", "RUNNING"):
            if row.get("execution_mode") != "executed" or not row.get("started_at"):
                raise ValueError("Executed event requires actual start timestamp")
            started = timestamp(row["started_at"])
            if status != "RUNNING":
                runtime = row.get("runtime_seconds")
                if (not row.get("finished_at") or timestamp(row["finished_at"]) < started
                        or isinstance(runtime, bool) or not isinstance(runtime, (float, int)) or not 0 <= runtime < float("inf")):
                    raise ValueError("Completed/failed execution requires ordered timestamps and measured runtime")
            if status == "COMPLETED" and not row.get("outputs_sha256"):
                raise ValueError("Completed execution requires an output digest")
            if status == "RUNNING" and (row.get("finished_at") is not None or row.get("runtime_seconds") is not None):
                raise ValueError("Running events cannot claim completed execution timing")
        else:
            if row.get("started_at") is not None or row.get("runtime_seconds") is not None:
                raise ValueError("Skipped/imported stages must not invent originating execution timings")
            if row.get("execution_mode") != ("imported" if status == "IMPORTED" else "skipped") or not row.get("reason") or not row.get("finished_at"):
                raise ValueError("Skipped/imported events require honest mode, observation timestamp and reason")
            timestamp(row["finished_at"])
        for field in ("input_artifacts", "output_artifacts"):
            for artifact in row.get(field, []):
                if not artifact.get("path") or not isinstance(artifact.get("size_bytes"), int) or artifact["size_bytes"] < 0 or len(artifact.get("sha256", "")) != 64:
                    raise ValueError("Artifact hashes and recorded byte sizes required")
    return value


class ExecutionJournal:
    """One journal per fresh campaign; repeated stages have distinct event IDs.

    `execute` runs the callable exactly once and re-raises failures after writing
    a sanitized failed event. Independent callers choose which stages continue.
    `skip` cannot emit COMPLETED or pretend to have execution timing.
    `import_artifacts` verifies bytes but explicitly does not claim model work.
    Callers must supply scientific input data only, never credentials/transports.
    """
    SKIP_STATES = {"SKIPPED_DEPENDENCY", "SKIPPED_UNCONFIGURED", "SKIPPED_NO_CANDIDATES", "SKIPPED_NO_INPUT"}
    EVIDENCE_STATES = {"NOT_ASSESSED", "AVAILABLE", "INSUFFICIENT_EVIDENCE", "FAILED_VALIDATION", "VALIDATED"}

    def __init__(self, path=None, *, run_id=None):
        self.path = Path(path).resolve() if path else None
        if self.path and self.path.exists():
            raise FileExistsError("A fresh execution journal must not overwrite an existing run")
        self.data = {"schema_version": "1.0", "run_id": run_id or str(uuid.uuid4()),
                     "started_at": utc_now(), "finished_at": None, "status": "RUNNING", "stages": []}
        self._timers = {}
        self._persist()

    @property
    def stages(self):
        return self.data["stages"]

    def _persist(self):
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = self.path.with_name(self.path.name + ".tmp")
            with temporary.open("w", encoding="utf8") as handle:
                json.dump(self.data, handle, sort_keys=True, indent=2, allow_nan=False)
                handle.flush()
                os.fsync(handle.fileno())
            temporary.replace(self.path)

    def _event(self, stage, *, inputs=None, dependencies=()):
        if self.data["status"] != "RUNNING":
            raise RuntimeError("Cannot append to a finalized execution journal")
        row = {"id": f"{len(self.stages)+1:04d}:{stage}", "stage": stage,
               "execution_status": "RUNNING", "evidence_status": "NOT_ASSESSED",
               "started_at": None, "finished_at": None, "runtime_seconds": None,
               "inputs_sha256": payload_hash(inputs), "outputs_sha256": None,
               "input_artifacts": [], "output_artifacts": [], "dependencies": list(dependencies),
               "reason": None, "error_type": None, "execution_mode": "executed"}
        self.stages.append(row)
        return row

    def start(self, stage, *, inputs=None, input_artifacts=(), dependencies=(), tools=()):
        row = self._event(stage, inputs=inputs, dependencies=dependencies)
        row["tools"] = json_value(list(tools))
        row["started_at"] = utc_now()
        self._timers[row["id"]] = time.perf_counter()
        self._persist()  # A process crash leaves a genuine RUNNING event.
        try:
            row["input_artifacts"] = hash_artifacts(input_artifacts)
            self._persist()
        except Exception as exc:
            self.fail(row, exc)
            raise
        return row

    def finish(self, event, *, outputs=None, output_artifacts=(), evidence_status="AVAILABLE", reason=None):
        if evidence_status not in self.EVIDENCE_STATES:
            raise ValueError("Unknown evidence status")
        if event["id"] not in self._timers or event["execution_status"] != "RUNNING":
            raise RuntimeError("Stage was not running")
        try:
            event["outputs_sha256"] = payload_hash(outputs)
            event["output_artifacts"] = hash_artifacts(output_artifacts)
        except Exception as exc:
            self.fail(event, exc)
            raise
        event.update(execution_status="COMPLETED", evidence_status=evidence_status, reason=reason,
                     finished_at=utc_now(), runtime_seconds=max(0.0, time.perf_counter()-self._timers.pop(event["id"])))
        self._persist()
        return event

    def fail(self, event, error):
        if event["id"] not in self._timers:
            return event
        event.update(execution_status="FAILED", evidence_status="INSUFFICIENT_EVIDENCE",
                     finished_at=utc_now(), runtime_seconds=max(0.0, time.perf_counter()-self._timers.pop(event["id"])),
                     error_type=type(error).__name__, reason="Stage execution failed; diagnostic exception type recorded without response bodies or secrets")
        self._persist()
        return event

    def execute(self, stage, function: Callable, *, inputs=None, input_artifacts=(), output_artifacts=(),
                dependencies=(), evidence_status="AVAILABLE", tools=()):
        row = self.start(stage, inputs=inputs, input_artifacts=input_artifacts, dependencies=dependencies, tools=tools)
        try:
            result = function()
            status = evidence_status(result) if callable(evidence_status) else evidence_status
            artifacts = output_artifacts(result) if callable(output_artifacts) else output_artifacts
            self.finish(row, outputs=result, output_artifacts=artifacts, evidence_status=status)
            return result
        except Exception as exc:
            self.fail(row, exc)
            raise

    def skip(self, stage, status, reason, *, inputs=None, dependencies=(), evidence_status="INSUFFICIENT_EVIDENCE"):
        if status not in self.SKIP_STATES or evidence_status not in self.EVIDENCE_STATES:
            raise ValueError("Skipped stages require an explicit skip and evidence status")
        row = self._event(stage, inputs=inputs, dependencies=dependencies)
        row.update(execution_status=status, evidence_status=evidence_status, execution_mode="skipped",
                   finished_at=utc_now(), reason=reason)
        self._persist()
        return row

    def import_artifacts(self, stage, paths, *, inputs=None, evidence_status="NOT_ASSESSED"):
        if evidence_status not in self.EVIDENCE_STATES:
            raise ValueError("Unknown evidence status")
        artifacts = hash_artifacts(paths)
        row = self._event(stage, inputs=inputs)
        row.update(execution_status="IMPORTED", evidence_status=evidence_status, execution_mode="imported",
                   finished_at=utc_now(), input_artifacts=artifacts, outputs_sha256=payload_hash(artifacts),
                   reason="Existing artifact bytes verified; no originating model execution claimed")
        self._persist()
        return row

    def latest(self, stage):
        return next((r for r in reversed(self.stages) if r["stage"] == stage), None)

    def finalize(self):
        if self._timers:
            raise RuntimeError("Cannot finalize while stages are running")
        if self.data["status"] != "RUNNING":
            return self.data
        failed = any(r["execution_status"] == "FAILED" for r in self.stages)
        skipped = any(r["execution_status"].startswith("SKIPPED") for r in self.stages)
        self.data.update(status="COMPLETED_WITH_FAILURES" if failed else "COMPLETED_WITH_SKIPS" if skipped else "COMPLETED",
                         finished_at=utc_now())
        self._persist()
        return self.data
