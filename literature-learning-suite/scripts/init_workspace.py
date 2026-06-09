#!/usr/bin/env python3
"""Initialize a portable literature workspace."""

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
DATA_ASSETS_DIR = SKILL_ROOT / "assets" / "data"
BUNDLED_DATA_FILES = (
    "bioc_genes.json",
    "kegg_pathways.json",
    "data-manifest.json",
)


def seed_workspace_data(root: Path, refresh: bool = False) -> list[str]:
    """Copy bundled runtime dictionaries into a workspace without clobbering edits."""
    target_dir = root / "data"
    target_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for name in BUNDLED_DATA_FILES:
        source = DATA_ASSETS_DIR / name
        target_name = "bundled-data-manifest.json" if name == "data-manifest.json" else name
        target = target_dir / target_name
        if not source.exists():
            raise FileNotFoundError(f"Missing bundled data asset: {source}")
        if refresh or not target.exists() or target.stat().st_size == 0:
            shutil.copy2(source, target)
            copied.append(target_name)
    return copied


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("literature-workspace"))
    parser.add_argument(
        "--refresh-data",
        action="store_true",
        help="replace workspace copies of bundled gene/pathway dictionaries",
    )
    args = parser.parse_args()
    root = args.root.expanduser().resolve()

    for name in (
        "fulltext",
        "fulltext_cache",
        "reports",
        "daily_digest",
        "data",
        "biorxiv_api",
        "cache",
        "imports",
        "exports",
        "logs",
        "config",
    ):
        (root / name).mkdir(parents=True, exist_ok=True)
    for name in (
        "papers.db",
        "concepts.db",
        "edges.db",
        "queries.db",
        "journal_metrics.db",
    ):
        (root / name).touch(exist_ok=True)

    copied_data = seed_workspace_data(root, refresh=args.refresh_data)

    manifest = root / "workspace.json"
    manifest_data = {}
    if manifest.exists():
        manifest_data = json.loads(manifest.read_text(encoding="utf-8-sig"))
    manifest_data.update(
        {
            "format": "literature-learning-suite",
            "version": 2,
            "bundled_data_manifest": "data/bundled-data-manifest.json",
        }
    )
    manifest_data.setdefault("created_at", datetime.now(timezone.utc).isoformat())
    manifest_data["updated_at"] = datetime.now(timezone.utc).isoformat()
    manifest.write_text(
        json.dumps(manifest_data, indent=2) + "\n",
        encoding="utf-8",
    )

    monitor_template = SKILL_ROOT / "assets" / "templates" / "monitor-job.json"
    monitor_config = root / "config" / "monitor-job.json"
    if monitor_template.exists() and not monitor_config.exists():
        shutil.copyfile(monitor_template, monitor_config)
    print(json.dumps({"root": str(root), "seeded_data": copied_data}, indent=2))


if __name__ == "__main__":
    main()
