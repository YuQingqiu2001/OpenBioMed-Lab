# Bundled Runtime Data

This directory contains the small reference datasets required for the biomedical edge
generator to work immediately after cloning.

- `bioc_genes.json`: 90,125 unique human and mouse gene symbols.
- `kegg_pathways.json`: 25,939 unique KEGG pathway and GO Biological Process labels.
- `data-manifest.json`: record counts, SHA-256 checksums, provenance, and regeneration path.

`scripts/init_workspace.py` copies the two runtime dictionaries into
`<workspace>/data/`. Existing workspace copies are preserved unless `--refresh-data` is used.
`scripts/gen_edges.py` also falls back to these bundled files when workspace copies are absent.

Run:

```bash
python scripts/check_assets.py
Rscript scripts/export_bioc_genes.R ./literature-workspace
```

The R exporter requires the listed Bioconductor packages and can refresh the dictionaries from
their upstream annotation sources. Source-package attribution and usage terms still apply.
