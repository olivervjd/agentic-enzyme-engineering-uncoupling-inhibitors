"""Explicitly configured OpenAI Responses transport, with no credential persistence."""
from __future__ import annotations

import json
import os
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from .gpt_rosalind import GPTRosalindJSONBackend
from ..evals.judges import RUBRIC

DEFAULT_EVIDENCE_MODEL = "gpt-rosalind-research"


class UnavailableModelTransport:
    """Preserve an explicit unavailable model without substituting curated synthesis."""
    def __init__(self, model, reason):
        self.model, self.reason = model, reason

    def __call__(self, operation, payload):
        raise RuntimeError(self.reason)


def configured_evidence_transport(model, api_key, *, curated_only=False):
    if curated_only:
        return None
    if not model or not model.strip():
        raise ValueError("An evidence model is required unless --curated-evidence-only is explicitly selected")
    if not api_key:
        raise ValueError("Evidence model requested but API key unavailable; use --use-codex-api-key or OPENAI_API_KEY")
    return OpenAIJSONTransport(model, api_key=api_key)


def object_schema(properties):
    return {"type": "object", "properties": properties, "required": list(properties), "additionalProperties": False}


def response_schema(operation):
    strings = {"type": "array", "items": {"type": "string"}}
    if operation in {"workflow_review", "workflow_judge"}:
        return object_schema({"summary": {"type": "string"}, "issues": strings, "next_actions": strings,
                              "recommendation": {"type": "string", "enum": ["revise", "fail", "ready_for_review"]}})
    if operation == "candidate_review":
        return object_schema({"assumptions": strings, "uncertainty": strings,
                              "recommendation": {"type": "string", "enum": ["advance", "reject", "more_computation"]}})
    if operation == "evidence_synthesis":
        return object_schema({"claims": {"type": "array", "items": object_schema({
            "text": {"type": "string"}, "source_ids": strings})}})
    if operation != "domain_judge":
        raise ValueError("Unsupported reasoning operation")
    return object_schema({"overall_score": {"type": "number"},
                          "category_scores": object_schema({k: {"type": "integer", "minimum": 1, "maximum": 5} for k in RUBRIC}),
                          **{k: strings for k in ("major_issues", "minor_issues", "unsafe_or_overclaimed_statements", "missing_evidence")},
                          "recommendation": {"type": "string", "enum": ["pass", "revise", "fail"]}})


class OpenAIJSONTransport:
    def __init__(self, model, api_key=None, timeout=180, opener=urlopen):
        if not model or not model.strip():
            raise ValueError("An explicit API model ID is required; no automatic model substitution")
        self.model = model
        self._key = api_key or os.getenv("OPENAI_API_KEY")
        if not self._key:
            raise ValueError("API key unavailable; supply OPENAI_API_KEY or explicitly reuse a Codex API-key login")
        self.timeout, self.opener = timeout, opener

    def __call__(self, operation, payload):
        body = {"model": self.model, "store": False,
                "instructions": payload.get("system_prompt", "Return only evidence-grounded structured output."),
                "input": json.dumps({k: v for k, v in payload.items() if k != "system_prompt"}, allow_nan=False),
                "text": {"format": {"type": "json_schema", "name": operation, "strict": True,
                                      "schema": response_schema(operation)}}}
        request = Request("https://api.openai.com/v1/responses", data=json.dumps(body).encode(),
                          headers={"Authorization": "Bearer " + self._key, "Content-Type": "application/json"})
        try:
            with self.opener(request, timeout=self.timeout) as response:
                result = json.load(response)
        except HTTPError as exc:
            code = exc.code
            exc.close()
            raise RuntimeError(f"OpenAI {operation} HTTP {code}; verify model access, quota and credentials") from None
        if result.get("status") != "completed":
            raise RuntimeError("OpenAI response did not complete; no substitute review generated")
        blocks = [b for item in result.get("output", []) if item.get("type") == "message" for b in item.get("content", [])]
        if any(b.get("type") == "refusal" for b in blocks):
            raise RuntimeError("Model declined the review")
        text = "".join(b.get("text", "") for b in blocks if b.get("type") == "output_text")
        value = json.loads(text)
        if not isinstance(value, dict):
            raise ValueError("Expected structured JSON object")
        return value


def configured_reasoner(model, api_key=None):
    backend = GPTRosalindJSONBackend(OpenAIJSONTransport(model, api_key=api_key))
    backend.model_id = model
    return backend
