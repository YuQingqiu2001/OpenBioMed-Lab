#!/usr/bin/env python3
"""Fetch legally accessible scholarly full text and cache normalized text."""

from __future__ import annotations

import argparse
import html.parser
import json
import re
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

from ll_common import http_get, safe_filename, utc_now


class HTMLTextExtractor(html.parser.HTMLParser):
    SKIP = {"script", "style", "nav", "header", "footer", "aside", "noscript"}

    def __init__(self):
        super().__init__()
        self.skip_depth = 0
        self.parts: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() in self.SKIP:
            self.skip_depth += 1

    def handle_endtag(self, tag):
        if tag.lower() in self.SKIP and self.skip_depth:
            self.skip_depth -= 1
        elif tag.lower() in {"p", "div", "section", "article", "h1", "h2", "h3", "li"}:
            self.parts.append("\n")

    def handle_data(self, data):
        if not self.skip_depth and data.strip():
            self.parts.append(data.strip())

    def value(self):
        return re.sub(r"\n{3,}", "\n\n", " ".join(self.parts)).strip()


def xml_text(raw: bytes) -> str:
    root = ET.fromstring(raw)
    chunks = []
    for tag in ("article-title", "abstract", "body"):
        for node in root.findall(f".//{tag}"):
            value = " ".join("".join(node.itertext()).split())
            if value:
                chunks.append(value)
    return "\n\n".join(dict.fromkeys(chunks))


def pubmed_to_pmc(pmid: str) -> str:
    params = urllib.parse.urlencode(
        {"dbfrom": "pubmed", "db": "pmc", "id": pmid, "retmode": "xml"}
    )
    raw, _ = http_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?" + params
    )
    root = ET.fromstring(raw)
    node = root.find(".//LinkSetDb/Link/Id")
    return node.text.strip() if node is not None and node.text else ""


def fetch_pubmed(pmid: str) -> tuple[str, str]:
    pmcid = pubmed_to_pmc(pmid)
    if not pmcid:
        raise RuntimeError("No PubMed Central full text is linked to this PMID")
    params = urllib.parse.urlencode(
        {"db": "pmc", "id": pmcid, "retmode": "xml"}
    )
    raw, _ = http_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    )
    return xml_text(raw), f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmcid}/"


def fetch_html(url: str) -> tuple[str, str]:
    raw, content_type = http_get(url, accept="text/html,application/xhtml+xml")
    if "xml" in content_type:
        return xml_text(raw), url
    parser = HTMLTextExtractor()
    parser.feed(raw.decode("utf-8", errors="replace"))
    value = parser.value()
    if len(value) < 500:
        raise RuntimeError("Fetched page did not contain enough readable text")
    return value, url


def save(root: Path, record_id: str, content: str, source_url: str) -> Path:
    folder = root / "fulltext"
    folder.mkdir(parents=True, exist_ok=True)
    stem = safe_filename(record_id)
    target = folder / f"{stem}.txt"
    target.write_text(content, encoding="utf-8")
    metadata = {
        "id": record_id,
        "source_url": source_url,
        "retrieved_at": utc_now(),
        "characters": len(content),
    }
    (folder / f"{stem}.metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return target


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("literature-workspace"))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pmid")
    group.add_argument("--arxiv")
    group.add_argument("--url")
    parser.add_argument("--id", dest="record_id")
    args = parser.parse_args()
    if args.pmid:
        content, source_url = fetch_pubmed(args.pmid)
        record_id = args.record_id or f"PMID:{args.pmid}"
    elif args.arxiv:
        source_url = f"https://arxiv.org/html/{args.arxiv}"
        content, source_url = fetch_html(source_url)
        record_id = args.record_id or f"ARXIV:{args.arxiv}"
    else:
        content, source_url = fetch_html(args.url)
        record_id = args.record_id or args.url
    print(save(args.root.expanduser().resolve(), record_id, content, source_url))


if __name__ == "__main__":
    main()
