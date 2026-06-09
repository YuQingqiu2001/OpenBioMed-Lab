#!/usr/bin/env python3
"""Export verified paper records as BibTeX, RIS, or Markdown."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from ll_common import normalize_doi, read_ndjson


def key(paper):
    author = paper.get("authors", ["unknown"])[0] if paper.get("authors") else "unknown"
    family = re.sub(r"[^A-Za-z0-9]", "", author.split()[-1]).lower() or "unknown"
    year = str(paper.get("year") or "0000")
    word = re.sub(r"[^A-Za-z0-9]", "", paper.get("title", "").split()[0]).lower()
    return f"{family}_{year}_{word or 'paper'}"


def bibtex(paper):
    fields = [
        f"  title = {{{paper.get('title', '')}}}",
        f"  author = {{{' and '.join(paper.get('authors', []))}}}",
        f"  year = {{{paper.get('year') or ''}}}",
    ]
    if paper.get("venue"):
        fields.append(f"  journal = {{{paper['venue']}}}")
    if normalize_doi(paper.get("doi")):
        fields.append(f"  doi = {{{normalize_doi(paper['doi'])}}}")
    if paper.get("url"):
        fields.append(f"  url = {{{paper['url']}}}")
    entry_type = "misc" if paper.get("source") == "arxiv" else "article"
    return f"@{entry_type}{{{key(paper)},\n" + ",\n".join(fields) + "\n}"


def ris(paper):
    lines = ["TY  - JOUR", f"TI  - {paper.get('title', '')}"]
    lines.extend(f"AU  - {author}" for author in paper.get("authors", []))
    lines.append(f"PY  - {paper.get('year') or ''}")
    if paper.get("venue"):
        lines.append(f"JO  - {paper['venue']}")
    if normalize_doi(paper.get("doi")):
        lines.append(f"DO  - {normalize_doi(paper['doi'])}")
    if paper.get("url"):
        lines.append(f"UR  - {paper['url']}")
    lines.append("ER  -")
    return "\n".join(lines)


def markdown(paper):
    authors = ", ".join(paper.get("authors", []))
    return (
        f"- {authors}. **{paper.get('title', '')}**. "
        f"{paper.get('venue', '')} ({paper.get('year') or ''}). "
        f"{normalize_doi(paper.get('doi')) or paper.get('id', '')}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", choices=("bibtex", "ris", "markdown"), default="bibtex")
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()
    papers, errors = read_ndjson(args.input)
    if errors:
        raise SystemExit(f"Input contains {len(errors)} invalid NDJSON lines")
    formatter = {"bibtex": bibtex, "ris": ris, "markdown": markdown}[args.format]
    separator = "\n\n" if args.format != "markdown" else "\n"
    content = separator.join(formatter(paper) for paper in papers) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
        print(args.output)
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
