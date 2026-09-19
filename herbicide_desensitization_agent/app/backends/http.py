from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Protocol
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


class JsonTransport(Protocol):
    def request(self, method: str, path: str, payload: dict[str, Any] | None = None) -> Any: ...


@dataclass
class UrllibJsonTransport:
    base_url: str
    api_key: str | None = None
    timeout_seconds: float = 120.0

    def __post_init__(self) -> None:
        if not self.base_url.startswith(("http://", "https://")):
            raise ValueError("NIM base URL must use HTTP or HTTPS")
        self.base_url = self.base_url.rstrip("/")

    def request(self, method: str, path: str, payload: dict[str, Any] | None = None) -> Any:
        body = json.dumps(payload).encode() if payload is not None else None
        headers = {"Accept": "application/json"}
        if body is not None:
            headers["Content-Type"] = "application/json"
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        request = Request(self.base_url + path, data=body, headers=headers, method=method.upper())
        try:
            with urlopen(request, timeout=self.timeout_seconds) as response:
                content = response.read()
        except HTTPError as exc:
            detail = exc.read().decode(errors="replace")[:1000]
            raise RuntimeError(f"NIM request failed with HTTP {exc.code}: {detail}") from exc
        except URLError as exc:
            raise RuntimeError(f"NIM request failed: {exc.reason}") from exc
        return json.loads(content) if content else {}

