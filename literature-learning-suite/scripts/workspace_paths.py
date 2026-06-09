"""Shared path resolution for the Hermes-compatible runtime."""

import os
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parent.parent
KG_ROOT = Path(
    os.environ.get("LITERATURE_KG_ROOT", str(Path.cwd() / "literature-workspace"))
).expanduser().resolve()
JOURNAL_METRICS_PATH = Path(
    os.environ.get(
        "JOURNAL_METRICS_PATH",
        str(KG_ROOT / "journal_metrics.db"),
    )
).expanduser().resolve()
DATA_ASSETS_DIR = SKILL_ROOT / "assets" / "data"
BUNDLED_JOURNAL_METRICS = DATA_ASSETS_DIR / "journal_metrics_2024.json"
BUNDLED_DATA_MANIFEST = DATA_ASSETS_DIR / "data-manifest.json"


def ensure_workspace() -> None:
    for path in (
        KG_ROOT,
        KG_ROOT / "data",
        KG_ROOT / "fulltext",
        KG_ROOT / "fulltext_cache",
        KG_ROOT / "reports",
        KG_ROOT / "daily_digest",
        KG_ROOT / "biorxiv_api",
    ):
        path.mkdir(parents=True, exist_ok=True)
