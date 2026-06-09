#!/usr/bin/env python
"""
Download bioRxiv/medRxiv metadata and abstracts via the official API.

This does NOT bypass Cloudflare to fetch PDF/JATS full text. It reliably downloads
API records: title, authors, DOI, date, category, abstract, JATS URL field, etc.

Examples:
  python scripts/download_biorxiv_api.py --from-date 2026-06-01 --to-date 2026-06-01
  python scripts/download_biorxiv_api.py --doi 10.64898/2026.05.31.727600
  python scripts/download_biorxiv_api.py --server medrxiv --from-date 2026-06-01
  python scripts/download_biorxiv_api.py --from-date 2026-06-01 --max-pages 1 --out-dir ./tmp/biorxiv_test
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import date, datetime
from pathlib import Path
from typing import Any

API_ROOT = "https://api.biorxiv.org/details"
from workspace_paths import KG_ROOT

DEFAULT_OUT_DIR = KG_ROOT / "biorxiv_api"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Download bioRxiv/medRxiv API records")
    p.add_argument("--server", choices=["biorxiv", "medrxiv"], default="biorxiv")
    p.add_argument("--from-date", default=None, help="YYYY-MM-DD")
    p.add_argument("--to-date", default=None, help="YYYY-MM-DD; default: same as --from-date")
    p.add_argument("--doi", default=None, help="Single bioRxiv/medRxiv DOI, e.g. 10.64898/2026.05.31.727600")
    p.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    p.add_argument("--max-pages", type=int, default=None, help="Limit API pages for testing")
    p.add_argument("--sleep", type=float, default=0.2, help="Seconds between pages")
    p.add_argument("--proxy", default=None, help="Optional proxy, e.g. http://127.0.0.1:7890")
    return p.parse_args()


def validate_date(s: str) -> str:
    datetime.strptime(s, "%Y-%m-%d")
    return s


def safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")


def build_opener(proxy: str | None) -> urllib.request.OpenerDirector:
    handlers = []
    if proxy:
        handlers.append(urllib.request.ProxyHandler({"http": proxy, "https": proxy}))
    opener = urllib.request.build_opener(*handlers)
    opener.addheaders = [
        ("User-Agent", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"),
        ("Accept", "application/json,text/plain,*/*"),
    ]
    return opener


def fetch_json(opener: urllib.request.OpenerDirector, url: str, timeout: int = 30) -> dict[str, Any]:
    with opener.open(url, timeout=timeout) as r:
        raw = r.read().decode("utf-8", errors="replace")
    return json.loads(raw)


def normalize_record(rec: dict[str, Any]) -> dict[str, Any]:
    return {
        "server": rec.get("server", ""),
        "date": rec.get("date", ""),
        "title": rec.get("title", ""),
        "authors": rec.get("authors", ""),
        "category": rec.get("category", ""),
        "doi": rec.get("doi", ""),
        "version": rec.get("version", ""),
        "type": rec.get("type", ""),
        "license": rec.get("license", ""),
        "jatsxml": rec.get("jatsxml", ""),
        "published": rec.get("published", ""),
        "abstract": rec.get("abstract", ""),
        "author_corresponding": rec.get("author_corresponding", ""),
        "author_corresponding_institution": rec.get("author_corresponding_institution", ""),
        "funder": rec.get("funder", ""),
    }


def write_outputs(records: list[dict[str, Any]], out_dir: Path, server: str, start: str, end: str) -> tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = safe_filename(f"{server}_{start}_to_{end}")
    jsonl_path = out_dir / f"{stem}.jsonl"
    csv_path = out_dir / f"{stem}.csv"
    md_path = out_dir / f"{stem}.md"

    fields = [
        "server", "date", "title", "authors", "category", "doi", "version", "type",
        "license", "jatsxml", "published", "abstract", "author_corresponding",
        "author_corresponding_institution", "funder",
    ]
    with jsonl_path.open("w", encoding="utf-8", newline="") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow(r)

    with md_path.open("w", encoding="utf-8") as f:
        f.write(f"# {server} records: {start} to {end}\n\n")
        f.write(f"Total downloaded: {len(records)}\n\n")
        for i, r in enumerate(records, 1):
            f.write(f"## {i}. {r.get('title','').strip()}\n\n")
            f.write(f"- DOI: {r.get('doi','')}\n")
            f.write(f"- Date: {r.get('date','')}\n")
            f.write(f"- Category: {r.get('category','')}\n")
            f.write(f"- Authors: {r.get('authors','')}\n")
            f.write(f"- JATS URL field: {r.get('jatsxml','')}\n\n")
            abs_text = (r.get("abstract") or "").replace("\r", " ").replace("\n", " ").strip()
            f.write(abs_text + "\n\n")
    return jsonl_path, csv_path, md_path


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    opener = build_opener(args.proxy)

    if args.doi:
        doi = args.doi.strip()
        url = f"{API_ROOT}/{args.server}/{doi}/na/json"
        try:
            payload = fetch_json(opener, url)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"ERROR fetching {url}: {e}", file=sys.stderr)
            return 2
        msg = (payload.get("messages") or [{}])[0]
        if str(msg.get("status", "")).lower() != "ok":
            print(f"API status not ok: {msg}", file=sys.stderr)
            return 3
        records = [normalize_record(x) for x in payload.get("collection", [])]
        print(f"doi={doi} records={len(records)}")
        jsonl_path, csv_path, md_path = write_outputs(records, out_dir, args.server, f"doi_{doi}", f"doi_{doi}")
        print("WROTE")
        print(f"jsonl={jsonl_path}")
        print(f"csv={csv_path}")
        print(f"md={md_path}")
        return 0

    if not args.from_date:
        print("ERROR: provide --from-date for date-range download, or --doi for a single paper", file=sys.stderr)
        return 1
    start = validate_date(args.from_date)
    end = validate_date(args.to_date or args.from_date)

    records: list[dict[str, Any]] = []
    cursor = 0
    page = 0
    total = None

    while True:
        page += 1
        url = f"{API_ROOT}/{args.server}/{start}/{end}/{cursor}/json"
        try:
            payload = fetch_json(opener, url)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"ERROR fetching {url}: {e}", file=sys.stderr)
            return 2

        msg = (payload.get("messages") or [{}])[0]
        if str(msg.get("status", "")).lower() != "ok":
            print(f"API status not ok: {msg}", file=sys.stderr)
            return 3

        batch = [normalize_record(x) for x in payload.get("collection", [])]
        records.extend(batch)
        total = int(msg.get("total") or len(records))
        print(f"page={page} cursor={cursor} batch={len(batch)} downloaded={len(records)} total={total}")

        if not batch:
            break
        if len(records) >= total:
            break
        if args.max_pages is not None and page >= args.max_pages:
            break
        cursor += len(batch)
        time.sleep(args.sleep)

    jsonl_path, csv_path, md_path = write_outputs(records, out_dir, args.server, start, end)
    print("WROTE")
    print(f"jsonl={jsonl_path}")
    print(f"csv={csv_path}")
    print(f"md={md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
