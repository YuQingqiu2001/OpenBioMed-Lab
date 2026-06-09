"""Shared utilities for Literature Learning Suite scripts."""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


USER_AGENT = "literature-learning-suite/1.2 (+https://github.com/)"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def http_get(url: str, accept: str | None = None, timeout: int = 45) -> tuple[bytes, str]:
    headers = {"User-Agent": USER_AGENT}
    if accept:
        headers["Accept"] = accept
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read(), response.headers.get("Content-Type", "")


def http_json(url: str, timeout: int = 45) -> dict[str, Any]:
    raw, _ = http_get(url, accept="application/json", timeout=timeout)
    return json.loads(raw)


def normalize_doi(value: str | None) -> str:
    doi = (value or "").strip().lower()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
    doi = re.sub(r"^doi:\s*", "", doi)
    return doi.rstrip(".,;")


def normalize_title(value: str | None) -> str:
    title = (value or "").casefold()
    title = re.sub(r"<[^>]+>", " ", title)
    title = re.sub(r"[^\w]+", " ", title, flags=re.UNICODE)
    return " ".join(title.split())


def stable_id(record: dict[str, Any]) -> str:
    if record.get("id"):
        return str(record["id"]).strip()
    if record.get("pmid"):
        return f"PMID:{record['pmid']}".strip()
    if record.get("arxiv_id"):
        return f"ARXIV:{record['arxiv_id']}".strip()
    doi = normalize_doi(record.get("doi"))
    if doi:
        return f"DOI:{doi}"
    fingerprint = normalize_title(record.get("title")) + "|" + str(record.get("year", ""))
    if fingerprint.strip("|"):
        return "HASH:" + hashlib.sha256(fingerprint.encode("utf-8")).hexdigest()[:16]
    return ""


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned[:180] or "record"


def read_ndjson(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    if not path.exists():
        return records, errors
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError("record is not a JSON object")
                records.append(value)
            except (json.JSONDecodeError, ValueError) as exc:
                errors.append({"line": line_number, "error": str(exc)})
    return records, errors


def append_ndjson(path: Path, records: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")
            count += 1
    return count


def atomic_write_ndjson(path: Path, records: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="\n", delete=False, dir=path.parent
    ) as handle:
        temp_path = Path(handle.name)
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")
    temp_path.replace(path)


def load_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8-sig"))
