#!/usr/bin/env python3
"""Import and query a user-supplied 2024 journal metrics dataset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from ll_common import atomic_write_ndjson, normalize_title, read_ndjson, utc_now


FIELD_ALIASES = {
    "journal": "name",
    "journal_name": "name",
    "title": "name",
    "abbreviation": "abbr_name",
    "impact_factor": "jif",
    "impact_factor_5y": "jif_5y",
    "five_year_impact_factor": "jif_5y",
    "category": "category_name",
    "jcr_quartile": "quartile",
}


def iter_source(path: Path) -> Iterable[dict[str, Any]]:
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as handle:
            yield from csv.DictReader(handle)
        return
    if suffix in {".jsonl", ".ndjson"}:
        records, errors = read_ndjson(path)
        if errors:
            raise ValueError(f"{path} contains {len(errors)} invalid NDJSON lines")
        yield from records
        return
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if isinstance(payload, list):
        yield from payload
    elif isinstance(payload, dict):
        records = payload.get("records") or payload.get("journals")
        if not isinstance(records, list):
            raise ValueError("JSON object must contain a records or journals array")
        yield from records
    else:
        raise ValueError("Expected a JSON array, JSON object, CSV, or NDJSON file")


def as_number(value: Any) -> float | None:
    if value in (None, "", "N/A", "n/a", "-", "—"):
        return None
    try:
        return float(str(value).replace(",", "").strip())
    except ValueError:
        return None


def normalize_record(record: dict[str, Any], metric_year: int) -> dict[str, Any]:
    renamed = {}
    for key, value in record.items():
        key_text = str(key).strip().casefold()
        normalized_key = FIELD_ALIASES.get(key_text, key_text)
        renamed[normalized_key] = value
    name = str(renamed.get("name", "")).strip()
    if not name:
        raise ValueError("journal name is required")
    output = {
        "name": name,
        "abbr_name": str(renamed.get("abbr_name", "")).strip(),
        "issn": str(renamed.get("issn", "")).strip(),
        "eissn": str(renamed.get("eissn", "")).strip(),
        "jif": as_number(renamed.get("jif")),
        "jif_5y": as_number(renamed.get("jif_5y")),
        "category_name": str(renamed.get("category_name", "")).strip(),
        "quartile": str(renamed.get("quartile", "")).strip(),
        "rank": renamed.get("rank"),
        "rank_total": renamed.get("rank_total"),
        "metric_year": metric_year,
    }
    return output


def import_metrics(source: Path, root: Path, metric_year: int, source_name: str) -> dict:
    records = []
    invalid = 0
    for raw in iter_source(source):
        try:
            records.append(normalize_record(raw, metric_year))
        except (TypeError, ValueError):
            invalid += 1
    if not records:
        raise ValueError("No valid journal metric records were found")
    root.mkdir(parents=True, exist_ok=True)
    target = root / "journal_metrics.db"
    atomic_write_ndjson(target, records)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    metadata = {
        "metric_year": metric_year,
        "record_count": len(records),
        "invalid_count": invalid,
        "source_name": source_name,
        "source_file_sha256": digest,
        "imported_at": utc_now(),
        "redistribution_note": (
            "The dataset was supplied locally and is not distributed with this skill. "
            "The operator is responsible for data licensing and attribution."
        ),
    }
    (root / "journal_metrics.metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


def lookup(root: Path, query: str, limit: int) -> list[dict[str, Any]]:
    records, errors = read_ndjson(root / "journal_metrics.db")
    if errors:
        raise ValueError("journal_metrics.db contains invalid NDJSON")
    needle = normalize_title(query)
    exact = query.casefold().replace("-", "")
    matches = []
    for record in records:
        identifiers = {
            str(record.get("issn", "")).casefold().replace("-", ""),
            str(record.get("eissn", "")).casefold().replace("-", ""),
        }
        names = " ".join(
            [str(record.get("name", "")), str(record.get("abbr_name", ""))]
        )
        if exact in identifiers or needle in normalize_title(names):
            matches.append(record)
        if len(matches) >= limit:
            break
    return matches


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("literature-workspace"))
    sub = parser.add_subparsers(dest="command", required=True)
    importer = sub.add_parser("import")
    importer.add_argument("source", type=Path)
    importer.add_argument("--year", type=int, default=2024)
    importer.add_argument("--source-name", default="user-supplied journal metrics")
    finder = sub.add_parser("lookup")
    finder.add_argument("query")
    finder.add_argument("-n", "--limit", type=int, default=10)
    sub.add_parser("status")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    if args.command == "import":
        if args.year != 2024:
            raise SystemExit("This release profile is for 2024 metrics; pass --year 2024")
        result = import_metrics(
            args.source.expanduser().resolve(), root, args.year, args.source_name
        )
    elif args.command == "lookup":
        result = lookup(root, args.query, max(1, args.limit))
    else:
        metadata = root / "journal_metrics.metadata.json"
        result = (
            json.loads(metadata.read_text(encoding="utf-8"))
            if metadata.exists()
            else {"installed": False, "metric_year": 2024}
        )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
