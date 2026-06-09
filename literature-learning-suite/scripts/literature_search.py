#!/usr/bin/env python3
"""Search public scholarly APIs with Python's standard library."""

import argparse
import json
import os
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta, timezone


USER_AGENT = "literature-learning-suite/1.2"


def request(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.read()


def text(element, path, namespaces=None):
    found = element.find(path, namespaces or {})
    return "".join(found.itertext()).strip() if found is not None else ""


def search_pubmed(query: str, limit: int):
    email = os.environ.get("PUBMED_EMAIL", "")
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": limit,
        "retmode": "json",
        "sort": "relevance",
        "tool": "literature-learning-suite",
    }
    if email:
        params["email"] = email
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    data = json.loads(request(base + "esearch.fcgi?" + urllib.parse.urlencode(params)))
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    fetch = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "xml",
    }
    if email:
        fetch["email"] = email
    root = ET.fromstring(request(base + "efetch.fcgi?" + urllib.parse.urlencode(fetch)))
    records = []
    for article in root.findall(".//PubmedArticle"):
        pmid = text(article, ".//PMID")
        title = text(article, ".//ArticleTitle")
        abstract = " ".join(
            "".join(x.itertext()).strip() for x in article.findall(".//AbstractText")
        )
        authors = []
        for author in article.findall(".//Author"):
            name = " ".join(
                value for value in (text(author, "ForeName"), text(author, "LastName"))
                if value
            )
            if name:
                authors.append(name)
        doi = ""
        for item in article.findall(".//ArticleId"):
            if item.attrib.get("IdType") == "doi":
                doi = (item.text or "").strip()
        records.append(
            {
                "id": f"PMID:{pmid}",
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "authors": authors,
                "venue": text(article, ".//Journal/Title"),
                "abstract": abstract,
                "source": "pubmed",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )
    return records


def search_arxiv(query: str, limit: int):
    ns = {"a": "http://www.w3.org/2005/Atom"}
    params = {
        "search_query": query,
        "start": 0,
        "max_results": limit,
        "sortBy": "relevance",
        "sortOrder": "descending",
    }
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode(params)
    root = ET.fromstring(request(url))
    records = []
    for entry in root.findall("a:entry", ns):
        raw_id = text(entry, "a:id", ns).split("/abs/")[-1]
        records.append(
            {
                "id": f"ARXIV:{raw_id}",
                "arxiv_id": raw_id,
                "title": " ".join(text(entry, "a:title", ns).split()),
                "authors": [
                    text(author, "a:name", ns)
                    for author in entry.findall("a:author", ns)
                ],
                "abstract": " ".join(text(entry, "a:summary", ns).split()),
                "source": "arxiv",
                "url": f"https://arxiv.org/abs/{raw_id}",
            }
        )
    return records


def search_crossref(query: str, limit: int):
    email = os.environ.get("CROSSREF_EMAIL", "")
    params = {
        "query.bibliographic": query,
        "rows": limit,
        "select": "DOI,title,author,published,container-title,URL,type,abstract",
    }
    if email:
        params["mailto"] = email
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    payload = json.loads(request(url))
    records = []
    for item in payload.get("message", {}).get("items", []):
        doi = item.get("DOI", "")
        date_parts = item.get("published", {}).get("date-parts", [[]])
        year = date_parts[0][0] if date_parts and date_parts[0] else None
        authors = []
        for author in item.get("author", []):
            name = " ".join(
                value for value in (author.get("given", ""), author.get("family", ""))
                if value
            )
            if name:
                authors.append(name)
        records.append(
            {
                "id": f"DOI:{doi.lower()}",
                "doi": doi,
                "title": " ".join(item.get("title", [""])),
                "authors": authors,
                "venue": " ".join(item.get("container-title", [""])),
                "abstract": item.get("abstract", ""),
                "year": year,
                "document_type": item.get("type", ""),
                "source": "crossref",
                "url": item.get("URL") or f"https://doi.org/{doi}",
            }
        )
    return records


def search_preprint(
    server: str,
    query: str,
    limit: int,
    from_date: str | None = None,
    to_date: str | None = None,
):
    end = date.fromisoformat(to_date) if to_date else datetime.now(timezone.utc).date()
    start = date.fromisoformat(from_date) if from_date else end - timedelta(days=30)
    words = [word.casefold() for word in query.split() if word.strip()]
    records = []
    cursor = 0
    while len(records) < limit:
        interval = f"{start.isoformat()}/{end.isoformat()}"
        url = f"https://api.biorxiv.org/details/{server}/{interval}/{cursor}/json"
        payload = json.loads(request(url))
        collection = payload.get("collection", [])
        if not collection:
            break
        for item in collection:
            haystack = f"{item.get('title', '')} {item.get('abstract', '')}".casefold()
            if words and not all(word in haystack for word in words):
                continue
            doi = item.get("doi", "")
            records.append(
                {
                    "id": f"{server.upper()}:{doi}",
                    "doi": doi,
                    "title": item.get("title", ""),
                    "authors": [
                        value.strip()
                        for value in item.get("authors", "").split(";")
                        if value.strip()
                    ],
                    "abstract": item.get("abstract", ""),
                    "source": server,
                    "url": f"https://www.{server}.org/content/{doi}",
                    "published": item.get("date", ""),
                    "version": item.get("version", ""),
                    "category": item.get("category", ""),
                }
            )
            if len(records) >= limit:
                break
        cursor += len(collection)
        total = int(payload.get("messages", [{}])[0].get("total", cursor))
        if cursor >= total:
            break
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source", choices=("pubmed", "arxiv", "crossref", "biorxiv", "medrxiv")
    )
    parser.add_argument("query")
    parser.add_argument("-n", "--limit", type=int, default=10)
    parser.add_argument("-o", "--output")
    parser.add_argument("--from-date", help="YYYY-MM-DD; preprint sources only")
    parser.add_argument("--to-date", help="YYYY-MM-DD; preprint sources only")
    args = parser.parse_args()
    limit = max(1, min(args.limit, 100))

    if args.source == "pubmed":
        records = search_pubmed(args.query, limit)
    elif args.source == "arxiv":
        records = search_arxiv(args.query, limit)
    elif args.source == "crossref":
        records = search_crossref(args.query, limit)
    else:
        records = search_preprint(
            args.source, args.query, limit, args.from_date, args.to_date
        )

    timestamp = datetime.now(timezone.utc).isoformat()
    for record in records:
        record["retrieved_at"] = timestamp
        record["provenance"] = {"query": args.query, "source": args.source}
    output = "\n".join(json.dumps(x, ensure_ascii=False) for x in records)
    if output:
        output += "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(output)
    else:
        sys.stdout.write(output)


if __name__ == "__main__":
    main()
