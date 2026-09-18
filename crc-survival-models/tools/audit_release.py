"""Create a reproducible structure, privacy and large-file audit for release."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".md", ".py", ".json", ".tsv", ".csv", ".txt", ".yml", ".yaml"}
REQUIRED = [
    "README.md",
    "SKILL.md",
    "docs/QUICKSTART.md",
    "docs/FIG3_PROVENANCE.md",
    "docs/THIRD_PARTY_NOTICES.md",
    "docs/RELEASE_CHECKLIST.md",
    "docs/validation/END_TO_END_VALIDATION.json",
    "LICENSE",
    "environment.yml",
    "predict.py",
    "predict_wsi.py",
    "models/cell_topology/MODEL_CARD.md",
    "models/cell_topology/model.py",
    "models/cell_topology/weights/cell_topology_model.npz",
    "models/patch_feature/MODEL_CARD.md",
    "models/patch_feature/model.py",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def main() -> None:
    files = sorted(path for path in ROOT.rglob("*") if path.is_file() and "__pycache__" not in path.parts)
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    for fold in range(5):
        for name in ("patch_encoder.pt", "spatial_field.pt", "spatial_transforms.pkl", "survival_aggregator.pt", "inference_contract.pkl", "training_manifest.json"):
            path = ROOT / f"models/patch_feature/weights/fold_{fold}/{name}"
            if not path.is_file():
                missing.append(relative(path))

    label_hits, absolute_hits, reference_identifier_headers = [], [], []
    banned = re.compile(r"(?i)(?:\bV5\b|\bV3\.3\b|v3_3)")
    absolute = re.compile(r"(?i)(?:[A-Z]:\\研究|/mnt/[a-z]/研究)")
    for path in files:
        if (
            path.suffix.lower() not in TEXT_SUFFIXES
            or "vendor" in path.parts
            or path == Path(__file__).resolve()
            or path.name == "PACKAGE_AUDIT.json"
        ):
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for number, line in enumerate(text.splitlines(), start=1):
            if banned.search(line):
                label_hits.append({"file": relative(path), "line": number, "text": line[:240]})
            if absolute.search(line):
                absolute_hits.append({"file": relative(path), "line": number, "text": line[:240]})
        if "reference_results" in path.parts and path.suffix.lower() in {".tsv", ".csv"}:
            header = text.splitlines()[0].lower() if text.splitlines() else ""
            if any(term in header.split("\t") for term in ("case_id", "sample_id", "patient_id")):
                reference_identifier_headers.append(relative(path))

    large = [{"file": relative(path), "bytes": path.stat().st_size} for path in files if path.stat().st_size >= 100_000_000]
    forbidden_third_party = [
        relative(path)
        for path in files
        if path.suffix.lower() == ".pth"
        or path.relative_to(ROOT).parts[0] in {"external", "vendor", "weights"}
    ]
    notices = (ROOT / "docs/THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    upstream_links_present = all(
        link in notices
        for link in (
            "https://github.com/TIO-IKIM/CellViT-plus-plus",
            "https://github.com/mahmoodlab/UNI",
        )
    )

    release_dir = ROOT / "release"
    release_dir.mkdir(parents=True, exist_ok=True)
    contents_path = release_dir / "PACKAGE_CONTENTS.tsv"
    with contents_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "bytes"])
        for path in files:
            if path.name not in {"PACKAGE_CONTENTS.tsv", "SHA256SUMS.tsv", "PACKAGE_AUDIT.json"}:
                writer.writerow([relative(path), path.stat().st_size])

    checksum_targets = [
        path
        for path in files
        if path.suffix.lower() in {".pth", ".pt", ".pkl", ".npz"}
        and "__pycache__" not in path.parts
    ]
    checksums_path = release_dir / "SHA256SUMS.tsv"
    with checksums_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["sha256", "bytes", "path"])
        for path in checksum_targets:
            writer.writerow([sha256(path), path.stat().st_size, relative(path)])

    structure_pass = (
        not missing
        and not label_hits
        and not absolute_hits
        and not reference_identifier_headers
        and not large
        and not forbidden_third_party
        and upstream_links_present
    )
    report = {
        "schema": "crc_survival_release_audit",
        "staging_pass": structure_pass,
        "public_release_ready": structure_pass,
        "public_release_blockers": [] if structure_pass else ["one or more automated release checks failed"],
        "checks": {
            "required_files_missing": missing,
            "historical_public_label_hits": label_hits,
            "local_absolute_path_hits": absolute_hits,
            "reference_files_with_identifier_columns": reference_identifier_headers,
            "large_files": large,
            "forbidden_third_party_files": forbidden_third_party,
            "official_upstream_links_present": upstream_links_present,
            "patch_folds_present": 5,
        },
        "artifacts": {
            "contents": contents_path.name,
            "checksums": checksums_path.name,
        },
    }
    (release_dir / "PACKAGE_AUDIT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if structure_pass else 2)


if __name__ == "__main__":
    main()
