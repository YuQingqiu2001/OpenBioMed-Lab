#!/usr/bin/env python3
"""Normalize heterogeneous paper records into the suite data model."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from ll_common import normalize_doi, normalize_title, stable_id, utc_now


def normalize_authors(value: Any) -> list[str]:
    if isinstance(value, list):
        result = []
        for item in value:
            if isinstance(item, dict):
                name = item.get("name") or " ".join(
                    str(item.get(key, "")) for key in ("given", "family")
                )
            else:
                name = str(item)
            name = " ".join(name.split())
            if name:
                result.append(name)
        return result
    if isinstance(value, str):
        delimiter = ";" if ";" in value else ","
        return [" ".join(x.split()) for x in value.split(delimiter) if x.strip()]
    return []


def normalize_record(raw: dict[str, Any], source_hint: str = "") -> dict[str, Any]:
    source = str(raw.get("source") or source_hint or "unknown").lower()
    title = raw.get("title") or raw.get("name") or ""
    if isinstance(title, list):
        title = title[0] if title else ""
    record = {
        "title": " ".join(str(title).split()),
        "authors": normalize_authors(raw.get("authors") or raw.get("author")),
        "year": raw.get("year"),
        "venue": raw.get("venue") or raw.get("journal") or raw.get("publisher") or "",
        "doi": normalize_doi(raw.get("doi") or raw.get("DOI")),
        "pmid": str(raw.get("pmid") or raw.get("PMID") or "").strip(),
        "arxiv_id": str(raw.get("arxiv_id") or raw.get("arxiv") or "").strip(),
        "source": source,
        "url": raw.get("url") or raw.get("URL") or "",
        "abstract": raw.get("abstract") or raw.get("summary") or "",
        "document_type": raw.get("document_type") or raw.get("type") or "article",
        "fulltext_status": raw.get("fulltext_status") or "metadata_only",
        "retrieved_at": raw.get("retrieved_at") or utc_now(),
        "provenance": raw.get("provenance") or {"raw": raw},
    }
    if not record["year"]:
        date_value = str(raw.get("published") or raw.get("date") or "")
        if len(date_value) >= 4 and date_value[:4].isdigit():
            record["year"] = int(date_value[:4])
    record["id"] = stable_id(record)
    record["normalized_title"] = normalize_title(record["title"])
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--source", default="")
    args = parser.parse_args()
    output = args.output.open("w", encoding="utf-8") if args.output else sys.stdout
    try:
        with args.input.open(encoding="utf-8-sig") as handle:
            for line in handle:
                if not line.strip():
                    continue
                raw = json.loads(line)
                record = normalize_record(raw, args.source)
                output.write(json.dumps(record, ensure_ascii=False) + "\n")
    finally:
        if args.output:
            output.close()


if __name__ == "__main__":
    main()
