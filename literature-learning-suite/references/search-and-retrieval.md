# Search and Retrieval

## Query Design

Build concept blocks before source syntax:

```text
(primary concept OR synonym OR controlled term)
AND (mechanism OR context)
AND (optional study design)
NOT (explicit exclusion)
```

Keep a query ledger with source, query, date, filters, cursor/page, and returned count.

## PubMed

Use MeSH plus title/abstract terms. Useful tags:

- `[MeSH Terms]`
- `[Title/Abstract]`
- `[Author]`
- `[Journal]`
- `[Publication Type]`
- `[PDAT]` or `[EDAT]`

For monitoring, EDAT is often better for newly indexed records. For systematic reviews,
document the full query exactly and avoid adding filters that are not validated.

NCBI E-utilities:

- `esearch.fcgi`: find PMIDs.
- `efetch.fcgi`: retrieve XML metadata and abstracts.
- `elink.fcgi`: map PMID to PMC where available.

Set a contact email and respect NCBI rate limits. Use an API key for higher permitted rates.

## arXiv

Useful fields:

- `all:`
- `ti:`
- `au:`
- `abs:`
- `cat:`

Use `AND`, `OR`, and `ANDNOT`. Preserve the version suffix of the paper actually analyzed.
Space requests by at least three seconds unless the client handles rate limiting.

## bioRxiv and medRxiv

The public details API supports date windows and DOI lookups:

```text
https://api.biorxiv.org/details/biorxiv/YYYY-MM-DD/YYYY-MM-DD/cursor/json
https://api.biorxiv.org/details/medrxiv/YYYY-MM-DD/YYYY-MM-DD/cursor/json
```

The APIs provide metadata and abstracts, not a guarantee of full-text access.

## Crossref

Use Crossref for DOI metadata and bibliographic verification:

```text
https://api.crossref.org/works/{doi}
https://api.crossref.org/works?query.bibliographic={query}
```

Include a descriptive `User-Agent` and contact address when operating at scale.

## Selection

Rank relevance separately from evidence quality. A high-impact or highly cited paper may be
useful for context but still have weak causal evidence.
