import io
import json
import unittest
from unittest.mock import Mock
from urllib.error import HTTPError

from herbicide_desensitization_agent.app.backends.openai_json import (
    DEFAULT_EVIDENCE_MODEL, configured_evidence_transport,
)
from herbicide_desensitization_agent.examples.check_model_access import inspect_model_access


class EvidenceConfigurationTests(unittest.TestCase):
    def test_documented_rosalind_id_is_default(self):
        self.assertEqual(DEFAULT_EVIDENCE_MODEL, "gpt-rosalind-research")
        self.assertEqual(configured_evidence_transport(DEFAULT_EVIDENCE_MODEL, "test-only").model,
                         "gpt-rosalind-research")

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
