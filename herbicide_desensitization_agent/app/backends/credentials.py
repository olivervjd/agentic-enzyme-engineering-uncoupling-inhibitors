"""Opt-in reuse of an actual Codex API key; never reuse ChatGPT OAuth tokens."""
from __future__ import annotations

import json
import os
from pathlib import Path


def resolve_api_key(*, use_codex=False, codex_home=None):
    if os.getenv("OPENAI_API_KEY"):
        return os.environ["OPENAI_API_KEY"], "environment"
    if not use_codex:
        return None, "unavailable"
    home = Path(codex_home or os.getenv("CODEX_HOME", Path.home() / ".codex"))
    try:
        data = json.loads((home / "auth.json").read_text())
    except FileNotFoundError:
        return None, "codex_key_unavailable"
    except (OSError, ValueError):
        raise ValueError("Cannot read Codex authentication configuration") from None
    key = data.get("OPENAI_API_KEY")
    if data.get("auth_mode") in {None, "apikey"} and isinstance(key, str) and key.strip():
        return key, "codex_api_key"
    return None, "codex_oauth_not_api_key"
