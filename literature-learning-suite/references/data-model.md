# Data Model

All databases use one JSON object per UTF-8 line.

## papers.db

Required fields:

- `id`
- `title`
- `source`
- `retrieved_at`

Recommended fields:

- `authors`, `year`, `venue`, `doi`, `pmid`, `arxiv_id`
- `abstract`, `url`, `document_type`
- `fulltext_status`, `fulltext_source`
- `provenance`
- `analysis`
- `entities`

## concepts.db

```json
{
  "id": "CONCEPT:example",
  "name": "Example",
  "type": "mechanism",
  "definition": "",
  "source_papers": ["PMID:123"],
  "created_at": ""
}
```

## edges.db

```json
{
  "source": "PMID:123",
  "target": "DOI:10.x/example",
  "relation": "supports",
  "description": "Specific reason for this relation.",
  "provenance": "analyst",
  "created_at": ""
}
```

## queries.db

```json
{
  "source": "pubmed",
  "query": "",
  "executed_at": "",
  "result_count": 0,
  "parameters": {}
}
```

## journal_metrics.db

Optional normalized records imported from a locally supplied authorized dataset. This release
profile requires `metric_year: 2024`. See `journal-metrics-2024.md` and
`assets/schemas/journal-metric-2024.schema.json`.

## Hermes Compatibility

The canonical public format stores deep analysis under `analysis` and biological entities under
`entities`. The bundled Hermes-derived tools also read the validated legacy layout:

- `tier2_subquestions`
- `tier3_ces_chains`
- `tier4_mechanism_cascade`
- `tier5_hidden_axis`
- `tier6_concept_innovation`
- `tier7_cross_refs`
- top-level `genes`, `pathways`, `cell_types`, `diseases`, and `technologies`

Do not duplicate both layouts in newly created records. Write the canonical nested format; use
the compatibility readers to consume existing Hermes databases without migration.

## Schemas and Templates

Canonical JSON Schemas:

- `assets/schemas/paper.schema.json`
- `assets/schemas/concept.schema.json`
- `assets/schemas/edge.schema.json`
- `assets/schemas/query.schema.json`
- `assets/schemas/monitor-job.schema.json`
- `assets/schemas/workspace.schema.json`

Editable record templates live in `assets/templates/`. The paper template contains the complete
T1-T7 layout, evidence grade, biological entities, full-text provenance, and cross-paper
relations.

## Bundled Runtime Dictionaries

`scripts/init_workspace.py` seeds:

- `data/bioc_genes.json`
- `data/kegg_pathways.json`
- `data/bundled-data-manifest.json`

The manifest records source packages, record counts, and SHA-256 values. Workspace copies are
never overwritten unless initialization is called with `--refresh-data`.

## Mutation Rules

- Append new records atomically.
- Reject duplicate stable IDs by default.
- Preserve raw metadata under `provenance.raw`.
- Store analysis revisions as dated revisions when history matters.
- Never treat `.db` as a binary or SQLite file.
