"""Check model access and optionally review a calibration report using the actual API."""
import argparse
import json
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import quote

from ..app.backends.calibration import write, sha
from ..app.backends.credentials import resolve_api_key
from ..app.backends.openai_json import DEFAULT_EVIDENCE_MODEL, DEFAULT_JUDGE_MODEL, OpenAIJSONTransport


def inspect_model_access(model, key, opener=urlopen):
    """Read model metadata only. This does not establish successful inference."""
    request = Request("https://api.openai.com/v1/models/" + quote(model, safe=""),
                      headers={"Authorization": "Bearer " + key})
    try:
        with opener(request, timeout=30) as response:
            result = json.load(response)
        return {"metadata_accessible": result.get("id") == model, "inference_verified": False}
    except HTTPError as exc:
        code = exc.code
        exc.close()
        return {"metadata_accessible": False, "http_status": code, "inference_verified": False,
                "status": "ACCESS_UNAVAILABLE" if code in {403, 404} else "CHECK_FAILED"}


def inspect_model_inference(model, key, opener=urlopen):
    """Run one bounded non-scientific request; return redacted access metadata only."""
    body = {"model": model, "input": "Reply with OK only.", "store": False,
            "reasoning": {"effort": "none"}, "max_output_tokens": 16}
    request = Request("https://api.openai.com/v1/responses", data=json.dumps(body).encode(),
                      headers={"Authorization": "Bearer " + key, "Content-Type": "application/json"})
    try:
        with opener(request, timeout=45) as response:
            result = json.load(response)
        return {"inference_verified": result.get("status") == "completed",
                "status": result.get("status"), "response_model": result.get("model"),
                "output_tokens": result.get("usage", {}).get("output_tokens")}
    except HTTPError as exc:
        code = exc.code
        exc.close()
        return {"inference_verified": False, "status": "HTTP_ERROR", "http_status": code}
    except URLError:
        return {"inference_verified": False, "status": "NETWORK_UNAVAILABLE"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--use-codex-api-key", action="store_true")
    parser.add_argument("--judge-model", default=DEFAULT_JUDGE_MODEL)
    parser.add_argument("--evidence-model", default=DEFAULT_EVIDENCE_MODEL)
    parser.add_argument("--binding-report", type=Path)
    parser.add_argument("--smoke-inference", action="store_true",
                        help="Make one bounded non-scientific request per distinct configured model")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    key, source = resolve_api_key(use_codex=args.use_codex_api_key)
    if not key:
        raise ValueError("No API key available")
    request = Request("https://api.openai.com/v1/models", headers={"Authorization": "Bearer " + key})
    with urlopen(request, timeout=30) as response:
        models = {m["id"] for m in json.load(response)["data"]}
    report = {"credential_source": source, "judge": {"model": args.judge_model, "listed": args.judge_model in models},
              "evidence": {"model": args.evidence_model, "listed": args.evidence_model in models,
                           **inspect_model_access(args.evidence_model, key)},
              "scope": "Model access preflight only; workflow not run" if not args.binding_report else
                       "Standalone calibration review; not a candidate review or completed workflow judge stage"}
    if args.smoke_inference:
        checks = {model: inspect_model_inference(model, key)
                  for model in dict.fromkeys((args.evidence_model, args.judge_model))}
        report["evidence"].update(checks[args.evidence_model])
        report["judge"].update(checks[args.judge_model])
        report["scope"] += "; bounded non-scientific inference smoke test"
    if args.binding_report and args.judge_model in models:
        binding = json.loads(args.binding_report.read_text())
        report["binding_report_sha256"] = sha(args.binding_report)
        report["calibration_review"] = OpenAIJSONTransport(args.judge_model, api_key=key)("candidate_review", {
            "system_prompt": "Review this binding calibration report, not a mutation. Do not infer experimental success. "
                             "Missing matched experimental labels requires more_computation. Identify assumptions and uncertainty. "
                             "This review cannot override deterministic scientific gates.", "calibration": binding})
        report["judge"]["request_verified"] = True
    write(args.output, report)
    print(json.dumps({k: report[k] for k in ("credential_source", "judge", "evidence")}))


if __name__ == "__main__":
    main()
