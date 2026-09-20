import io
import json
import unittest
from unittest.mock import Mock, patch
from urllib.error import HTTPError

from herbicide_desensitization_agent.app.backends.openai_json import (
    DEFAULT_AGENT_MODEL, DEFAULT_EVIDENCE_MODEL, DEFAULT_REVIEW_MODEL, DEFAULT_JUDGE_MODEL,
    configured_evidence_transport, configured_reasoner, response_schema,
)
from herbicide_desensitization_agent.examples.check_model_access import inspect_model_access, inspect_model_inference


class EvidenceConfigurationTests(unittest.TestCase):
    def test_requested_luna_id_is_default(self):
        self.assertEqual(DEFAULT_EVIDENCE_MODEL, "gpt-5.6-luna")
        self.assertEqual(configured_evidence_transport(DEFAULT_EVIDENCE_MODEL, "test-only").model,
                         "gpt-5.6-luna")

    def test_no_implicit_curated_fallback_without_credentials(self):
        with self.assertRaisesRegex(ValueError, "API key unavailable"):
            configured_evidence_transport(DEFAULT_EVIDENCE_MODEL, None)

    def test_curated_mode_must_be_explicit(self):
        self.assertIsNone(configured_evidence_transport(None, None, curated_only=True))
        with self.assertRaisesRegex(ValueError, "evidence model"):
            configured_evidence_transport(None, "test-only")

    def test_model_access_rejection_is_not_enabled(self):
        response = io.BytesIO(b'{"error": {"message": "do not persist this body"}}')
        opener = Mock(side_effect=HTTPError("https://api.openai.com", 404, "Not found", {}, response))
        result = inspect_model_access(DEFAULT_EVIDENCE_MODEL, "test-only", opener)
        self.assertEqual(result["status"], "ACCESS_UNAVAILABLE")
        self.assertFalse(result["metadata_accessible"])
        self.assertFalse(result["inference_verified"])
        self.assertTrue(response.closed)
        self.assertNotIn("test-only", json.dumps(result))

    def test_metadata_success_is_not_inference_success(self):
        opener = Mock(return_value=io.BytesIO(json.dumps({"id": DEFAULT_EVIDENCE_MODEL}).encode()))
        result = inspect_model_access(DEFAULT_EVIDENCE_MODEL, "test-only", opener)
        self.assertTrue(result["metadata_accessible"])
        self.assertFalse(result["inference_verified"])
        self.assertEqual(opener.call_args.args[0].get_method(), "GET")

    def test_other_errors_not_mislabeled_missing_access(self):
        opener = Mock(side_effect=HTTPError("https://api.openai.com", 429, "Rate limit", {}, io.BytesIO()))
        self.assertEqual(inspect_model_access(DEFAULT_EVIDENCE_MODEL, "test-only", opener)["status"], "CHECK_FAILED")


class LunaMigrationTests(unittest.TestCase):
    @patch.dict("os.environ", {}, clear=True)
    def test_all_live_role_defaults_use_requested_model(self):
        from herbicide_desensitization_agent.examples.run_integrated_epsps import build_parser
        args = build_parser().parse_args(["--previous", "prior", "--calibration", "cal", "--output", "out", "--boltz-cache", "cache"])
        self.assertEqual({args.evidence_model, args.review_model, args.judge_model,
                          DEFAULT_AGENT_MODEL, DEFAULT_REVIEW_MODEL, DEFAULT_JUDGE_MODEL}, {"gpt-5.6-luna"})

    def test_legacy_imports_resolve_to_neutral_implementation(self):
        from herbicide_desensitization_agent.app.backends.gpt_rosalind import GPTRosalindJSONBackend
        from herbicide_desensitization_agent.app.backends.gpt_json import GPTJSONBackend
        from herbicide_desensitization_agent.app.backends.interfaces import ReasoningBackend, RosalindReasoningBackend
        self.assertIs(GPTRosalindJSONBackend, GPTJSONBackend)
        self.assertIs(RosalindReasoningBackend, ReasoningBackend)
        backend = configured_reasoner(DEFAULT_AGENT_MODEL, api_key="test-only")
        self.assertIsInstance(backend, GPTJSONBackend)
        self.assertEqual(backend.model_id, DEFAULT_AGENT_MODEL)

    def test_reasoning_outputs_cannot_add_measurement_fields(self):
        for operation in ("candidate_review", "evidence_synthesis", "workflow_review", "workflow_judge"):
            schema = response_schema(operation)
            self.assertFalse(schema["additionalProperties"])
            self.assertFalse({"affinity", "kd", "ki", "ddg", "rmsd", "decision"} & set(schema["properties"]))
        claims = response_schema("evidence_synthesis")["properties"]["claims"]["items"]
        self.assertEqual(set(claims["properties"]), {"text", "source_ids"})
        self.assertFalse(claims["additionalProperties"])

    def test_smoke_is_bounded_and_returns_no_content_or_key(self):
        opener = Mock(return_value=io.BytesIO(json.dumps({"status": "completed", "model": DEFAULT_AGENT_MODEL,
                         "output": [{"private": "not-exported"}], "usage": {"output_tokens": 1}}).encode()))
        report = inspect_model_inference(DEFAULT_AGENT_MODEL, "test-only", opener)
        self.assertTrue(report["inference_verified"])
        body = json.loads(opener.call_args.args[0].data)
        self.assertEqual(body["model"], DEFAULT_AGENT_MODEL)
        self.assertEqual(body["max_output_tokens"], 16)
        self.assertFalse(body["store"])
        for secret in ("test-only", "not-exported"):
            self.assertNotIn(secret, json.dumps(report))

    def test_smoke_failure_does_not_become_success_or_leak_body(self):
        response = io.BytesIO(b'{"error":"private-error-body"}')
        opener = Mock(side_effect=HTTPError("https://api.openai.com", 403, "Forbidden", {}, response))
        report = inspect_model_inference(DEFAULT_AGENT_MODEL, "test-only", opener)
        self.assertFalse(report["inference_verified"])
        self.assertEqual(report["http_status"], 403)
        self.assertTrue(response.closed)
        self.assertNotIn("private-error-body", json.dumps(report))
