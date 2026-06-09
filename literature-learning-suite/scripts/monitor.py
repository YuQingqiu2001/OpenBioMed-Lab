#!/usr/bin/env python3
"""Run a resumable multi-source monitoring job from JSON configuration."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ll_common import append_ndjson, read_ndjson, stable_id, utc_now
from literature_search import search_arxiv, search_crossref, search_preprint, search_pubmed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--root", type=Path, default=Path("literature-workspace"))
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    for name in ("papers.db", "queries.db"):
        (root / name).touch(exist_ok=True)
    existing, errors = read_ndjson(root / "papers.db")
    if errors:
        raise SystemExit("papers.db contains invalid NDJSON")
    ids = {stable_id(record) for record in existing}
    found = duplicate = added = failed = 0
    query_records = []
    new_records = []
    limit = int(config.get("batch_size", 10))
    for item in config.get("queries", []):
        if item.get("enabled", True) is False:
            continue
        source = item["source"]
        query = item["query"]
        started = utc_now()
        try:
            if source == "pubmed":
                records = search_pubmed(query, limit)
            elif source == "arxiv":
                records = search_arxiv(query, limit)
            elif source == "crossref":
                records = search_crossref(query, limit)
            elif source in ("biorxiv", "medrxiv"):
                records = search_preprint(
                    source,
                    query,
                    limit,
                    item.get("from_date"),
                    item.get("to_date"),
                )
            else:
                raise ValueError(f"Unsupported source: {source}")
            found += len(records)
            for record in records:
                record["retrieved_at"] = utc_now()
                record["provenance"] = {"source": source, "query": query}
                record_id = stable_id(record)
                record["id"] = record_id
                if record_id in ids:
                    duplicate += 1
                    continue
                ids.add(record_id)
                new_records.append(record)
                added += 1
            query_records.append(
                {
                    "source": source,
                    "query": query,
                    "executed_at": started,
                    "result_count": len(records),
                    "status": "success",
                }
            )
        except Exception as exc:
            failed += 1
            query_records.append(
                {
                    "source": source,
                    "query": query,
                    "executed_at": started,
                    "result_count": 0,
                    "status": "error",
                    "error": str(exc),
                }
            )
    append_ndjson(root / "papers.db", new_records)
    append_ndjson(root / "queries.db", query_records)
    manifest = {
        "completed_at": utc_now(),
        "found": found,
        "duplicates": duplicate,
        "added": added,
        "failed_queries": failed,
    }
    reports = root / "reports"
    reports.mkdir(exist_ok=True)
    target = reports / ("run-" + manifest["completed_at"].replace(":", "-") + ".json")
    target.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**manifest, "manifest": str(target)}, indent=2))
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
