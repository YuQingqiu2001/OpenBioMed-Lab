# Literature Learning Suite — Complete User Guide

> **Version**: 1.3.0 | **License**: CC BY-NC-SA 4.0

---

## Table of Contents

1. [Overview & Design Philosophy](#1-overview--design-philosophy)
2. [System Architecture](#2-system-architecture)
3. [Quick Start](#3-quick-start)
4. [Core Methodology: S-tier 7-Layer Protocol](#4-core-methodology-s-tier-7-layer-protocol)
5. [Knowledge Graph Data Model](#5-knowledge-graph-data-model)
6. [Edge Generation Engine v3.1](#6-edge-generation-engine-v31)
7. [Toolchain Reference](#7-toolchain-reference)
8. [Automated Monitoring](#8-automated-monitoring)
9. [Trap Library & Troubleshooting](#9-trap-library--troubleshooting)
10. [Platform Notes](#10-platform-notes)
11. [Best Practices](#11-best-practices)
12. [FAQ](#12-faq)
13. [Glossary](#13-glossary)

---

## 1. Overview & Design Philosophy

### 1.1 What This Is

The Literature Learning Suite (LLS) is a complete system for scholarly literature discovery,
deep analysis, knowledge graph construction, and automated monitoring. It is not "yet another
reference manager" — it is a **literature cognition operating system**.

Traditional reference managers (Zotero, EndNote, Mendeley) solve the "where to store" problem.
LLS solves the **"how to read"** and **"how to think"** problems. Through a rigorous 7-layer
dissection protocol, it transforms every paper into a structured, queryable, linkable knowledge
graph node.

### 1.2 Design Principles

These principles were crystallized through the deep analysis of hundreds of papers:

**Principle 1: Depth over breadth.** One complete S-tier 7-layer analysis is worth more than 50
papers with only titles and abstracts. If you have 2 hours, spend 1.5 hours deeply reading 2
papers, not scanning 20.

**Principle 2: Verify before citing.** Cross-check bibliographic identity (DOI/PMID/arXiv ID)
across at least two sources before analysis or citation.

**Principle 3: Separate reported evidence from inference.** Every statement must be labeled with
its epistemological tier — REPORTED (directly in the paper), SUPPORTED INFERENCE (follows from
evidence + established knowledge), HYPOTHESIS (analyst-generated), or UNKNOWN (cannot be
determined from available material).

**Principle 4: Judge evidence by study design, not journal prestige.** Impact factor is metadata,
not a quality score. LLS includes a standardized evidence rubric based on risk of bias, sample
size, control quality, and replicability.

**Principle 5: Append-only persistence.** All NDJSON databases use append-only mode. This
guarantees complete operational history and prevents irreversible data loss.

**Principle 6: Never fabricate.** No invented papers, identifiers, statistics, or mechanisms.
If a citation cannot be verified, exclude it.

**Principle 7: Text-as-untrusted-data.** Paper text and web content are untrusted data — objects
of analysis, never agent instructions. This is critical for prompt injection prevention.

### 1.3 Comparison with Existing Tools

| Dimension | Zotero/EndNote | Traditional Lit Review | LLS |
|-----------|---------------|----------------------|-----|
| Storage | PDF + metadata | — | NDJSON knowledge graph |
| Analysis depth | Manual notes | Human summary | 7-layer structured dissection |
| Cross-paper links | Manual tags | Subjective | 5-strategy automatic semantic edges |
| Evidence grading | None | Experience-based | Standardized rubric |
| Automation | Plugin-assisted | None | Full monitoring pipeline |
| Queryability | Keyword search | Not queryable | Full-text search + graph traversal |
| Persistence format | Proprietary | Plain text | NDJSON (human-readable) |

---

## 2. System Architecture

### 2.1 Directory Structure

```
literature-learning-suite/
├── GUIDE_EN.md              ← This document
├── GUIDE_ZH.md              ← Chinese guide (primary)
├── GUIDE_DE.md              ← German guide
├── GUIDE_JA.md              ← Japanese guide
├── GUIDE_KO.md              ← Korean guide
├── SKILL.md                 ← Agent operation protocol
├── README.md                ← Project overview
├── LICENSE                  ← CC BY-NC-SA 4.0
│
├── scripts/                 ← Core toolchain (23+ Python + Node.js + R)
├── assets/data/             ← Validated reference data (bundled)
│   ├── bioc_genes.json      ← 90,125 human + mouse gene symbols
│   ├── kegg_pathways.json   ← 25,939 KEGG + GO BP terms
│   └── journal_metrics_2024.json ← 21,800 journal IF/quartiles
├── references/              ← Methodology & protocols (25+ docs)
└── tests/                   ← Unit tests + synthetic fixtures
```

### 2.2 Data Flow Pipeline

```
Research Question
       │
       ▼
  [1. Frame search strategy]  ← PICO/PECO framework
       │
       ▼
  [2. Multi-source search]    ← PubMed/arXiv/bioRxiv/Crossref
       │
       ▼
  [3. Verify identity]        ← Cross-check DOI/PMID
       │
       ▼
  [4. Acquire full text]      ← PMC XML / arXiv HTML / CDP extraction
       │
       ▼
  [5. S-tier 7-layer analysis] ← LLM deep reasoning (T1-T7)
       │
       ▼
  [6. Grade evidence]         ← Risk of bias / sample size / replicability
       │
       ▼
  [7. Persist knowledge graph] ← NDJSON append
       │
  ┌────┼────┐
  ▼    ▼    ▼
 [Edges] [Digest] [Network]
       │
       ▼
  [9. Quality selfcheck + Monitor]
```

### 2.3 Runtime Workspace

Created by `init_workspace.py`, a self-contained directory:

```
my-workspace/
├── workspace.json           ← Metadata
├── papers.db                ← NDJSON, one paper per line
├── concepts.db              ← NDJSON, one concept per line
├── edges.db                 ← NDJSON, one semantic edge per line
├── queries.db               ← Search log
├── data/                    ← Seeded gene/pathway dictionaries
├── fulltext/                ← Downloaded documents
├── daily_digest/            ← Auto-generated summaries
├── config/                  ← Monitor configuration
├── exports/                 ← BibTeX exports
└── logs/                    ← Runtime logs
```

---

## 3. Quick Start

### 3.1 Requirements

- Python 3.10+
- pip
- Node.js v24+ (bioRxiv CDP extraction only)
- R + Bioconductor (gene dictionary regeneration only)

### 3.2 Installation

```bash
git clone <repo-url>
cd literature-learning-suite
pip install -r scripts/requirements.txt
```

### 3.3 Initialize Workspace

```bash
python scripts/init_workspace.py --root ./my-workspace
```

This creates the workspace directory structure and seeds the gene dictionary (90,125 symbols)
and pathway labels (25,939 terms).

### 3.4 Set Environment Variable

```bash
# Linux/macOS
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"

# Windows PowerShell
$env:LITERATURE_KG_ROOT = (Resolve-Path ./my-workspace)
```

If unset, scripts default to `./literature-workspace`.

### 3.5 First Search

```bash
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10
python scripts/literature_search.py arxiv "single cell foundation model" -n 10
python scripts/literature_search.py biorxiv "spatial transcriptomics" -n 10
```

### 3.6 Ingest Papers

```bash
python scripts/normalize_records.py search-results.jsonl -o normalized.jsonl
python scripts/validate_records.py normalized.jsonl
python scripts/kg.py --root ./my-workspace add normalized.jsonl
python scripts/kg.py --root ./my-workspace stats
```

---

## 4. Core Methodology: S-tier 7-Layer Protocol

This is the heart of LLS. Every research paper undergoes complete 7-layer analysis.

### 4.1 Why 7 Layers

Traditional reading stops at "read abstract → tag → write two sentences." This shallow
processing means you cannot capture deep logical structure, discover implicit connections
between papers, evaluate true evidence strength, or form queryable structured knowledge.

The 7-layer protocol decomposes reading into seven independent but complementary cognitive
dimensions, forcing the analyst (human or LLM) to produce substantive output in each.

### 4.2 Layer 1 (T1): Bibliographic & Study Profile

Establish verifiable identity and context: full title, first/corresponding authors and
affiliations, journal (auto-annotated with JCR 2024 IF and quartile via `enrich_paper_if()`),
stable identifiers, document type, study design, sample size, full-text status.

This layer can be automated.

### 4.3 Layer 2 (T2): Core Scientific Question

Distill the paper into **one falsifiable core question** and decompose into **≥5 testable
subquestions**. The core question must be mechanistic, not descriptive.

**Example (good):**
> Core: How does GDF15 induce NK cell dysfunction through xenobiotic receptor signaling?
> Subquestions: Which tumors express GDF15? Which receptors mediate the signal? What
> signaling cascade impairs NK killing? Is the dysfunction reversible? Does blocking
> GDF15 restore anti-tumor immunity in vivo?

### 4.4 Layer 3 (T3): Claim-Evidence-Synthesis Chains (≥5)

Transform core findings into ≥5 independent C-E-S chains with concrete evidence.

| Field | Requirement | Threshold |
|-------|------------|-----------|
| `claim` | Falsifiable conclusion in your own words | — |
| `evidence` | Concrete data: N, effect size, p-value, assay | **>20 chars** (empty-shell trigger) |
| `synthesis` | How evidence supports/weakens claim; cross-domain links | — |
| `strength` | 1-5 stars, justified by evidence quality | — |
| `uncertain` | Alternative explanations, missing validation | — |

**Good:** "2847 NSCLC patients, SHAP interaction analysis, EFS HR=0.52 (95% CI 0.34-0.78),
p=0.002, independent radiological review."

**Bad (empty shell):** "The authors proved their hypothesis."

### 4.5 Layer 4 (T4): Mechanism Cascade

Map the complete causal chain from trigger to phenotype, precise to modification sites:

```
Trigger → Receptor → Second messenger/kinase → Transcription factor → Target gene → Phenotype
```

Must include: ≥3 cascade steps, key modification sites (e.g., "NF-κB p65 Ser536
phosphorylation" not "NF-κB activation"), downstream effects, ≥1 feedback loop.
Label each step: `demonstrated`, `supported`, `background`, or `hypothesis`.

### 4.6 Layer 5 (T5): Hidden Organizing Axes (≥3)

Uncover patterns the paper does **not** explicitly state — implicit assumptions,
spatial/temporal organization, selection effects, or deeper experimental logic.

Each axis: an `observation` (verifiable fact from the paper) + an `interpretation`
(the deeper pattern you discovered).

T5 tests your synthesis ability, not your ability to copy the Discussion section.

### 4.7 Layer 6 (T6): Conceptual Contribution

Identify: new concepts with operational definitions (≥1), prior views challenged or narrowed
(≥1), methodological breakthroughs (≥1), and boundary conditions (when does this NOT
generalize?).

### 4.8 Layer 7 (T7): Cross-Paper Relations (≥5)

Connect the paper to other verified records using substantive biological relations
(`supports`, `contradicts`, `extends`, `replicates`, `shared_mechanism`, etc.).
Each relation requires a 60-150 word description explaining WHY.

**Explicitly prohibited:** `same_journal`, `same_issue`, `same_author` (non-biological noise).

### 4.9 Empty-Shell S Detection

A record labeled `S` is an empty shell if ANY of:
- tier2_subquestions empty
- tier3_ces_chains < 5
- Any evidence ≤ 20 chars
- tier4 cascade < 3 steps
- tier5_hidden_axis < 3
- tier7_cross_refs < 5
- Any tier7 relation is prohibited

**Post-write verification:**
```bash
python -c "
import json
with open('./my-workspace/papers.db', 'r', encoding='utf-8-sig') as f:
    for line in f:
        p = json.loads(line)
        if p.get('analysis_tier') != 'S': continue
        chains = p.get('tier3_ces_chains', [])
        empty_ev = sum(1 for c in chains if len(c.get('evidence','')) <= 20)
        print(f'{p[\"id\"]}: T2={len(p.get(\"tier2_subquestions\",[]))} '
              f'T3={len(chains)} empty_evidence={empty_ev}')
"
```

### 4.10 Non-Research Papers

Papers classified as news/editorial/erratum/retraction: label `analysis_tier: "NR"`,
record the actual type in `core_findings`, and do NOT fabricate analysis content.

---

## 5. Knowledge Graph Data Model

### 5.1 Design Rationale

LLS uses NDJSON (Newline-Delimited JSON) instead of SQLite because:
1. Human-readable and editable
2. Version-control-friendly (git diff works line-by-line)
3. Append-only for data integrity
4. Queryable with standard CLI tools (head, grep, jq)
5. Zero external database dependencies

### 5.2 papers.db

One complete paper record per line. Required fields: `id`, `title`, `source`, `retrieved_at`.
Full S-tier v4.0 schema includes tier2-tier7 analysis fields plus `entities` for
gene/pathway/cell_type/disease annotations.

### 5.3 concepts.db

```json
{"id": "CONCEPT:immunometabolic_checkpoint", "name": "Immunometabolic Checkpoint",
 "type": "mechanism", "definition": "...", "source_papers": ["PMID:42251595"]}
```

### 5.4 edges.db

```json
{"source": "PMID:42251595", "target": "PMID:39988000", "relation": "extends",
 "description": "60-150 word biological rationale", "provenance": "analyst"}
```

### 5.5 queries.db

Search audit log: `source`, `query`, `executed_at`, `result_count`, `parameters`.

---

## 6. Edge Generation Engine v3.1

### 6.1 Overview

`gen_edges.py` builds semantically meaningful connections between papers using 5 complementary
strategies, all implemented with inverted indices (O(M) complexity, not O(N²)).

### 6.2 Five Strategies

| # | Strategy | Relation | Logic |
|---|----------|----------|-------|
| 1 | Explicit refs | extends, accompanied_by | tier7_cross_refs from deep analysis (non-bio filtered) |
| 2 | Shared molecules | shares_molecules (≥2) | 90,125 genes, two-stage matching |
| 2.5 | Text overlap | shares_topic (≥4) | Core-findings bag-of-words → inverted index |
| 3 | Disease × method | shares_disease_method | Disease label × method label cross-product |
| 4 | Hidden axis | shares_paradigm | Tier5 deep-pattern keyword resonance |
| 5 | Concept nodes | defines_concept | Paper → concept edges from concepts.db |

### 6.3 Usage

```bash
# CRITICAL: clear bytecode cache before each run
rm -rf scripts/__pycache__
python -B scripts/gen_edges.py
python scripts/build_network.py  # refresh visualization
```

---

## 7. Toolchain Reference

### Search & Discovery

- `literature_search.py` — Multi-source search (pubmed/arxiv/crossref/biorxiv)
- `search_arxiv.py` — arXiv-specific with author/category filters
- `download_biorxiv_api.py` — bioRxiv/medRxiv batch download

### Verification & Normalization

- `verify_citation.py` — Cross-check DOI/PMID/arXiv identity
- `normalize_records.py` — Standardize + deduplicate
- `validate_records.py` — JSON Schema validation

### Full-Text Acquisition

- `fulltext_fetch.py` — Unified downloader (PMC XML + arXiv HTML)
- `extract_pymupdf.py` — PDF text + table extraction
- `extract_marker.py` — PDF OCR extraction
- `extract_biorxiv_cdp.mjs` — Chrome CDP bioRxiv/medRxiv extraction

### Knowledge Graph

- `kg.py` — CLI: add/stats/search/audit
- `kg_core.py` — Library: CRUD + automatic IF annotation via `enrich_paper_if()`

### Analysis & Output

- `gen_edges.py` — Semantic edge generator v3.1
- `gen_digest.py` — Daily markdown digest
- `build_network.py` — Interactive force-directed HTML graph
- `selfcheck_knowledge_graph.py` — 10-dimension quality audit
- `export_citations.py` — BibTeX/CSL export

### Automation

- `monitor.py` — Resumable batch monitor
- `init_workspace.py` — Bootstrap a new workspace

---

## 8. Automated Monitoring

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json --root ./my-workspace
```

Features: resumable, batch processing, idempotent (dedup by stable ID), fault-isolated.

For agent-host cron:
```json
{"name": "Daily Literature Deep Analysis", "schedule": "0 9 * * *",
 "skills": ["literature-learning-suite"], "enabled_toolsets": ["terminal", "file", "web"]}
```

⚠️ Avoid interactive tools in unattended cron runs. bioRxiv CDP extraction is manual.

---

## 9. Trap Library & Troubleshooting

### TRAP-001: .pyc Bytecode Ghost
**Symptom:** Source changes have no effect.
**Prevention:** `rm -rf scripts/__pycache__ && python -B`

### TRAP-002: Sandbox Write Isolation
**Symptom:** Files written in `execute_code` don't appear on disk.
**Prevention:** Use `write_file` tool or `terminal` + Python.

### TRAP-003: .db File Rejected by read_file
**Symptom:** `Cannot read binary file 'papers.db'`
**Prevention:** Use `terminal` + `python -c` or `head` to read .db files.

### TRAP-004: Unicode SyntaxError
**Symptom:** `SyntaxError: invalid character '→'`
**Prevention:** Replace `→` with `->`, `""` with `'`, `ΔΨ` with `Delta`/`Psi`.

### TRAP-005: arXiv Silent Failure
**Symptom:** Zero results, no error.
**Prevention:** Use `curl -sL "http://export.arxiv.org/api/query?..."` (HTTP + -L).

### TRAP-006: PubMed Proxy Slowdown
**Symptom:** PubMed queries time out.
**Prevention:** `unset http_proxy https_proxy` before PubMed calls.

### TRAP-007: Empty-Shell S Records
**Symptom:** `analysis_tier: "S"` but T2-T7 empty.
**Detection:** Run `selfcheck_knowledge_graph.py`, check `s_empty` count.

### TRAP-008: Non-Biological Edge Pollution
**Symptom:** `gen_edges.py` produces `same_journal` edges.
**Prevention:** Strategy 1 has built-in NON_BIO_RELS filter.

### TRAP-009: bioRxiv JATS XML 403
**Symptom:** `.source.xml` endpoint returns 403.
**Workaround:** Chrome CDP extraction or API abstracts (300-500 words).

### TRAP-010: Journal IF Fuzzy Match Error
**Symptom:** "Cell" IF matched to "Cell Reports".
**Fix:** Length ratio > 0.7 threshold in `kg_core.py`.

---

## 10. Platform Notes

### Linux / macOS
- `python3` is typical command
- Chrome CDP: `google-chrome --remote-debugging-port=9223`
- Gene regeneration: `Rscript scripts/export_bioc_genes.R`

### Windows
- Use git-bash (MSYS) or WSL
- `python` (not `python3`)
- Chrome CDP: `scripts/biorxiv_chrome_cdp_launcher.bat`
- Avoid non-ASCII in MSYS-executed Python one-liners (see TRAP-004)

### All Platforms
- Bundled gene/pathway dictionaries and JCR 2024 metrics work out of the box
- `enrich_paper_if()` auto-annotates IF from 21,800-journal database
- `.db` files are plain NDJSON text

---

## 11. Best Practices

### Search Strategy
- Start narrow, then broaden
- PubMed first (faster, ~780ms direct), then arXiv
- Log every search for reproducibility
- Zero results ≠ evidence of absence — check syntax, connectivity, filters

### Analysis Strategy
- Full text > abstract. Only-abstract → limit claims, label `abstract_only`
- No batching. One paper at a time. Template-copying is the #1 source of empty shells
- Evidence must be concrete: "2847 NSCLC, SHAP, EFS HR=0.52" not "authors proved"
- Mechanism precise to site: "NF-κB p65 Ser536 phosphorylation" not "activates pathway"

### Maintenance
- Daily `selfcheck_knowledge_graph.py` as last step of cron
- Clear .pyc after every gen_edges.py change
- Weekly backup of papers.db, edges.db, concepts.db
- Monitor logs/ for anomalies

---

## 12. FAQ

**Q: Why NDJSON instead of SQLite?**
A: Human-readable, git-friendly, CLI-queryable, zero dependencies. For "write once, read many,
occasionally append" workloads, NDJSON is a better fit.

**Q: How to handle incremental updates?**
A: `add_paper()` and `kg.py add` have built-in dedup by stable ID. Just append new records.

**Q: How often to run gen_edges.py?**
A: After every batch of new papers. In daily cron, run as last step.

**Q: Can I use my own journal IF data?**
A: Yes. Put NDJSON in `workspace/journal_metrics.db`. System prefers workspace data then falls
back to bundled JCR 2024 data.

**Q: How to get bioRxiv full text?**
A: API abstracts (300-500 words) are always available. JavaScript-rendered full text requires
Chrome CDP extraction with manual Cloudflare verification.

**Q: Chinese-language papers?**
A: Current search sources are English-focused. Chinese papers can be manually created and
ingested via `normalize_records.py` + `kg.py add`. The analysis protocol is language-agnostic.

**Q: How to upgrade legacy analysis records?**
A: See `references/s-tier-upgrade-workflow.md`. Flow: audit → batch upgrade (10/set) →
rebuild edges and network → re-audit.

---

## 13. Glossary

| Term | Definition |
|------|------------|
| Knowledge Graph | NDJSON-based network of papers, concepts, and semantic edges |
| S-tier Analysis | Complete 7-layer dissection protocol |
| Empty-shell S | Record labeled S with no substantive T2-T7 content |
| Claim-Evidence-Synthesis Chain (CES) | Core unit of T3 analysis |
| Hidden Organizing Axis | T5 layer — patterns the paper does not explicitly state |
| Semantic Edge | Biologically meaningful connection between papers |
| NDJSON | Newline-Delimited JSON storage format |
| Inverted Index | Core data structure of gen_edges.py (O(M) complexity) |
| Workspace | Runtime data directory (my-workspace/) |
| Stable Identifier | Persistently citable ID (PMID/DOI/arXiv ID) |
| Evidence Rubric | Standardized evidence grading system |
| CDP | Chrome DevTools Protocol |
| Two-stage Matching | Regex candidate extraction → set lookup confirmation |

---

> Maintained by the Literature Learning Suite project.
> Version 1.3.0, last updated 2026-06-09.
> License: CC BY-NC-SA 4.0
