#!/usr/bin/env python3
"""Verify DOI, PMID, or arXiv identifiers against public registries."""

import argparse
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


def fetch(url: str):
    req = urllib.request.Request(
        url, headers={"User-Agent": "literature-learning-suite/1.0"}
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.status, response.read()


def verify_doi(doi: str):
    status, raw = fetch("https://api.crossref.org/works/" + urllib.parse.quote(doi))
    item = json.loads(raw).get("message", {})
    return {
        "verified": status == 200,
        "source": "crossref",
        "doi": item.get("DOI", doi),
        "title": (item.get("title") or [""])[0],
        "authors": [
            " ".join(x for x in (a.get("given", ""), a.get("family", "")) if x)
            for a in item.get("author", [])
        ],
        "type": item.get("type", ""),
        "url": item.get("URL", ""),
    }


def verify_pmid(pmid: str):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    )
    status, raw = fetch(url)
    payload = json.loads(raw)
    item = payload.get("result", {}).get(pmid, {})
    return {
        "verified": status == 200 and bool(item),
        "source": "pubmed",
        "pmid": pmid,
        "title": item.get("title", ""),
        "authors": [a.get("name", "") for a in item.get("authors", [])],
        "venue": item.get("fulljournalname", ""),
        "doi": next(
            (
                x.get("value", "")
                for x in item.get("articleids", [])
                if x.get("idtype") == "doi"
            ),
            "",
        ),
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    }


def verify_arxiv(arxiv_id: str):
    url = (
        "https://export.arxiv.org/api/query?"
        + urllib.parse.urlencode({"id_list": arxiv_id})
    )
    status, raw = fetch(url)
    root = ET.fromstring(raw)
    ns = {"a": "http://www.w3.org/2005/Atom"}
    entry = root.find("a:entry", ns)
    return {
        "verified": status == 200 and entry is not None,
        "source": "arxiv",
        "arxiv_id": arxiv_id,
        "title": " ".join((entry.findtext("a:title", "", ns) if entry is not None else "").split()),
        "url": f"https://arxiv.org/abs/{arxiv_id}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--doi")
    group.add_argument("--pmid")
    group.add_argument("--arxiv")
    args = parser.parse_args()
    if args.doi:
        result = verify_doi(args.doi)
    elif args.pmid:
        result = verify_pmid(args.pmid)
    else:
        result = verify_arxiv(args.arxiv)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
