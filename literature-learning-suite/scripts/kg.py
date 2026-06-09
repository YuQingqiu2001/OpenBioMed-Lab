#!/usr/bin/env python3
"""Append, query, and audit an NDJSON literature library."""

import argparse
import json
from collections import Counter
from pathlib import Path

from ll_common import normalize_doi, stable_id


DATABASES = (
    "papers.db",
    "concepts.db",
    "edges.db",
    "queries.db",
    "journal_metrics.db",
)


def read_ndjson(path: Path):
    records = []
    errors = []
    if not path.exists():
        return records, errors
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                errors.append({"line": line_number, "error": str(exc)})
    return records, errors


def add(root: Path, input_path: Path):
    db = root / "papers.db"
    existing, errors = read_ndjson(db)
    if errors:
        raise SystemExit(f"Refusing append: papers.db has {len(errors)} invalid lines")
    ids = {stable_id(x) for x in existing if stable_id(x)}
    added = duplicate = invalid = 0
    with input_path.open(encoding="utf-8-sig") as source, db.open("a", encoding="utf-8") as out:
        for line in source:
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                invalid += 1
                continue
            record_id = stable_id(record)
            if not record_id or not record.get("title") or not record.get("source"):
                invalid += 1
                continue
            if record.get("doi"):
                record["doi"] = normalize_doi(record["doi"])
            record["id"] = record_id
            if record_id in ids:
                duplicate += 1
                continue
            out.write(json.dumps(record, ensure_ascii=False) + "\n")
            ids.add(record_id)
            added += 1
    print(json.dumps({"added": added, "duplicates": duplicate, "invalid": invalid}))


def stats(root: Path):
    result = {}
    for name in DATABASES:
        records, errors = read_ndjson(root / name)
        result[name] = {"records": len(records), "invalid_lines": len(errors)}
    papers, _ = read_ndjson(root / "papers.db")
    result["sources"] = Counter(x.get("source", "unknown") for x in papers)
    result["fulltext_status"] = Counter(
        x.get("fulltext_status", "unknown") for x in papers
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


def search(root: Path, query: str):
    words = [x.casefold() for x in query.split()]
    papers, _ = read_ndjson(root / "papers.db")
    for paper in papers:
        haystack = json.dumps(paper, ensure_ascii=False).casefold()
        if all(word in haystack for word in words):
            print(json.dumps(paper, ensure_ascii=False))


def audit(root: Path):
    report = {"databases": {}, "issues": []}
    for name in DATABASES:
        records, errors = read_ndjson(root / name)
        report["databases"][name] = {
            "records": len(records),
            "invalid_lines": errors,
        }
    papers, _ = read_ndjson(root / "papers.db")
    ids = [stable_id(x) for x in papers]
    duplicates = [key for key, count in Counter(ids).items() if key and count > 1]
    missing = [
        index + 1
        for index, paper in enumerate(papers)
        if not stable_id(paper) or not paper.get("title") or not paper.get("source")
    ]
    if duplicates:
        report["issues"].append({"duplicate_ids": duplicates})
    if missing:
        report["issues"].append({"missing_required_fields_at_records": missing})
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(1 if report["issues"] else 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("literature-workspace"))
    sub = parser.add_subparsers(dest="command", required=True)
    add_parser = sub.add_parser("add")
    add_parser.add_argument("input", type=Path)
    sub.add_parser("stats")
    search_parser = sub.add_parser("search")
    search_parser.add_argument("query")
    sub.add_parser("audit")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    for name in DATABASES:
        (root / name).touch(exist_ok=True)
    if args.command == "add":
        add(root, args.input)
    elif args.command == "stats":
        stats(root)
    elif args.command == "search":
        search(root, args.query)
    else:
        audit(root)


if __name__ == "__main__":
    main()
