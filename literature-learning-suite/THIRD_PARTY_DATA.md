# Third-Party Data Notice

The skill bundles generated lookup files so entity matching and journal annotation
work immediately after cloning:

| File | Contents | Records | Upstream sources |
|---|---:|---:|---|
| `assets/data/bioc_genes.json` | Human and mouse gene symbols | 90,125 | `org.Hs.eg.db`, `org.Mm.eg.db` |
| `assets/data/kegg_pathways.json` | Pathway and biological-process labels | 25,939 | `KEGGREST`, `GO.db` |
| `assets/data/journal_metrics_2024.json` | JCR journal impact factors | 21,800 | JCR 2024 (Clarivate) |

Exact SHA-256 values are stored in `assets/data/data-manifest.json`.

Gene Ontology content requires Gene Ontology attribution and follows its upstream
license.  Bioconductor annotation package terms continue to apply.  KEGG access and
redistribution terms may differ by use case; review the current KEGG terms before
redistributing a newly generated copy as part of another product.

JCR journal metrics are from Clarivate Journal Citation Reports 2024.  The bundled
file contains journal names, impact factors, quartiles, and category information
(identifiers and metrics), not article abstracts or full bibliographic records.
Review Clarivate redistribution terms for your use case before distributing
modified or newly generated copies in another product.

Regenerate local workspace copies with:

```bash
Rscript scripts/export_bioc_genes.R ./literature-workspace
```

The exporter uses upstream packages installed by the operator. It does not
download or bundle commercial journal metrics — those must be obtained from
Clarivate directly.
