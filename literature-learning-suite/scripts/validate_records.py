#!/usr/bin/env python3
"""Validate paper or edge NDJSON records without external dependencies."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ll_common import normalize_doi, read_ndjson, stable_id


ROOT = Path(__file__).resolve().parents[1]
ONTOLOGY_PATH = ROOT / "assets" / "data" / "relation-ontology.json"
ONTOLOGY = json.loads(ONTOLOGY_PATH.read_text(encoding="utf-8-sig"))
ALLOWED_RELATIONS = {item["id"] for item in ONTOLOGY["relations"]}
FORBIDDEN_RELATIONS = set(ONTOLOGY.get("forbidden_relations", []))
FULLTEXT_STATUSES = {"fulltext", "abstract_only", "metadata_only", "unavailable"}
ANALYSIS_TIERS = {"B", "A", "S"}


def validate_paper(record: dict) -> list[str]:
    issues = []
    for field in ("title", "source", "retrieved_at"):
        if not record.get(field):
            issues.append(f"missing {field}")
    if not stable_id(record):
        issues.append("missing stable identifier")
    if record.get("doi") and "/" not in normalize_doi(record["doi"]):
        issues.append("malformed DOI")
    if record.get("authors") is not None and not isinstance(record["authors"], list):
        issues.append("authors must be an array")
    if record.get("analysis") is not None and not isinstance(record["analysis"], dict):
        issues.append("analysis must be an object")
    if record.get("entities") is not None and not isinstance(record["entities"], dict):
        issues.append("entities must be an object")
    if (
        record.get("fulltext_status") is not None
        and record["fulltext_status"] not in FULLTEXT_STATUSES
    ):
        issues.append("invalid fulltext_status")
    if (
        record.get("analysis_tier") is not None
        and record["analysis_tier"] not in ANALYSIS_TIERS
    ):
        issues.append("analysis_tier must be B, A, or S")
    return issues


def validate_edge(record: dict) -> list[str]:
    issues = []
    for field in ("source", "target", "relation", "description"):
        if not record.get(field):
            issues.append(f"missing {field}")
    if record.get("source") == record.get("target"):
        issues.append("self-loop")
    relation = record.get("relation")
    if relation in FORBIDDEN_RELATIONS:
        issues.append("forbidden relation")
    elif relation and relation not in ALLOWED_RELATIONS:
        issues.append("unknown relation")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--type", choices=("paper", "edge"), default="paper")
    args = parser.parse_args()
    records, parse_errors = read_ndjson(args.input)
    failures = [{"line": x["line"], "issues": [x["error"]]} for x in parse_errors]
    invalid_records = 0
    validator = validate_paper if args.type == "paper" else validate_edge
    for index, record in enumerate(records, 1):
        issues = validator(record)
        if issues:
            invalid_records += 1
            failures.append({"record": index, "id": stable_id(record), "issues": issues})
    print(
        json.dumps(
            {
                "records": len(records),
                "parse_errors": len(parse_errors),
                "valid": len(records) - invalid_records,
                "failures": failures,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
