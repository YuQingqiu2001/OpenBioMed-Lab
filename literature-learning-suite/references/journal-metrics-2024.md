# Journal Metrics: 2024 Profile

This skill supports locally supplied journal metrics with `metric_year: 2024`.

## Publication Boundary

The public skill contains:

- A neutral JSON/CSV/NDJSON importer.
- A normalized schema for journal name, ISSN, JIF, five-year JIF, category, quartile, and rank.
- Lookup commands and provenance metadata.
- Small synthetic test fixtures.

The public skill does not contain a full Journal Citation Reports export or other licensed
commercial dataset. Public visibility of individual values does not by itself grant permission
to redistribute a compiled database.

Clarivate's public 2024 release page states that JCR covers more than 21,800 journals. Its
separate first-inclusion page offers a specific downloadable list of 544 journals, while access
to the full indicators is directed through the JCR product or an institutional subscription:

- https://clarivate.com/news/clarivate-reveals-worlds-leading-and-trusted-journals-with-the-2024-journal-citation-reports/
- https://clarivate.com/academia-government/first-time-journal-citation-reports-inclusion-list-2024/

## Import

```bash
python scripts/journal_metrics.py --root ./literature-workspace import \
  /path/to/authorized-journal-metrics-2024.json \
  --year 2024 \
  --source-name "Authorized 2024 source"
```

The import creates:

- `journal_metrics.db`: normalized UTF-8 NDJSON.
- `journal_metrics.metadata.json`: metric year, source label, SHA-256, record count, and import time.

## Lookup

```bash
python scripts/journal_metrics.py --root ./literature-workspace lookup "Nature"
python scripts/journal_metrics.py --root ./literature-workspace lookup "0028-0836"
python scripts/journal_metrics.py --root ./literature-workspace status
```

## Interpretation

Always call the value a **2024 journal metric** in outputs. Do not present journal-level metrics
as evidence that an individual paper is reliable, important, or methodologically sound.
