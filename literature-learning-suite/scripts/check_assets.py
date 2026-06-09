#!/usr/bin/env python3
"""Verify bundled JSON assets, runtime dictionaries, and declared checksums."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
DATA = ASSETS / "data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_assets() -> dict:
    json_files = sorted(ASSETS.rglob("*.json"))
    parsed = {}
    failures = []
    for path in json_files:
        try:
            parsed[path] = json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"{path.relative_to(ROOT)}: {exc}")

    manifest_path = DATA / "data-manifest.json"
    manifest = parsed.get(manifest_path, {})
    datasets = manifest.get("datasets", []) if isinstance(manifest, dict) else []
    verified = []
    for dataset in datasets:
        path = ROOT / dataset["path"]
        if not path.exists():
            failures.append(f"missing dataset: {dataset['path']}")
            continue
        value = parsed.get(path)
        if not isinstance(value, list):
            failures.append(f"{dataset['path']}: expected JSON array")
            continue
        actual_hash = sha256(path)
        actual_count = len(value)
        if actual_hash != dataset["sha256"]:
            failures.append(f"{dataset['path']}: checksum mismatch")
        if actual_count != dataset["record_count"]:
            failures.append(f"{dataset['path']}: record count mismatch")
        # Dedup check: serialize dict entries via JSON for hashing
        try:
            unique = len(set(json.dumps(v, sort_keys=True, ensure_ascii=False) for v in value))
        except TypeError:
            unique = len(set(str(v) for v in value))
        if unique != actual_count:
            failures.append(
                f"{dataset['path']}: duplicate values ({actual_count} total, {unique} unique)"
            )
        verified.append(
            {
                "path": dataset["path"],
                "records": actual_count,
                "sha256": actual_hash,
            }
        )

    return {
        "json_files": len(json_files),
        "datasets": verified,
        "failures": failures,
    }


def main() -> None:
    report = check_assets()
    print(json.dumps(report, indent=2))
    raise SystemExit(1 if report["failures"] else 0)


if __name__ == "__main__":
    main()
