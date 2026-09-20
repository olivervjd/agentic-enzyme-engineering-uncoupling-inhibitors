"""Bounded Europe PMC abstract retrieval with reproducible, untrusted source records."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen

from .calibration import write


class EuropePMCRetriever:
    def __init__(self, output, *, limit=8, opener=urlopen):
        if not 1 <= limit <= 25:
            raise ValueError("Literature limit must be between 1 and 25")
        self.root, self.limit, self.opener = Path(output), limit, opener

    def retrieve(self, entry):
        query = f'("{entry.protein_name}" OR "{entry.agi}") AND "{entry.herbicide}"'
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({
            "query": query, "format": "json", "resultType": "core", "pageSize": self.limit,
        })
        with self.opener(Request(url, headers={"Accept": "application/json"}), timeout=45) as response:
            raw = response.read(4_000_001)
        if len(raw) > 4_000_000:
            raise ValueError("Literature response exceeds size limit")
        data = json.loads(raw)
        records = data["resultList"]["result"]
        sources, seen = [], set()
        for record in records[:self.limit]:
            identity = str(record["source"]) + ":" + str(record["id"])
            if identity in seen or not record.get("abstractText"):
                continue
            seen.add(identity)
            source_id = "epmc-" + hashlib.sha256(identity.encode()).hexdigest()[:16]
            artifact = self.root / (source_id + ".json")
            write(artifact, record)
            sources.append({
                "id": source_id, "url": "https://europepmc.org/article/" +
                quote(str(record["source"]), safe="") + "/" + quote(str(record["id"]), safe=""),
                "artifact": str(artifact.resolve()), "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
                "title": record.get("title", ""), "abstract": record["abstractText"][:20000],
                "scope": "Retrieved abstract; organism, endpoint and target applicability require review",
                "evidence_level": "abstract_only_unreviewed",
            })
        self.root.mkdir(parents=True, exist_ok=True)
        (self.root / "search_response.json").write_bytes(raw)
        report = {"status": "COMPLETED" if sources else "NO_RESULTS", "query": query,
                  "url": url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                  "response_sha256": hashlib.sha256(raw).hexdigest(), "sources": sources,
                  "legend": "Bounded relevance-ranked Europe PMC search, abstracts only. Retrieval is not "
                            "experimental validation or an exhaustive review. No calibration labels are inferred."}
        write(self.root / "retrieval.json", report)
        return report
