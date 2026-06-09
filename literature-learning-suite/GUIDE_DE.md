# Literature Learning Suite — Vollständiges Benutzerhandbuch

> **Version**: 1.3.0 | **Lizenz**: CC BY-NC-SA 4.0

---

## Inhaltsverzeichnis

1. [Überblick & Designphilosophie](#1-überblick--designphilosophie)
2. [Systemarchitektur](#2-systemarchitektur)
3. [Schnellstart](#3-schnellstart)
4. [Kernmethodik: S-Tier-7-Schichten-Protokoll](#4-kernmethodik-s-tier-7-schichten-protokoll)
5. [Wissensgraph-Datenmodell](#5-wissensgraph-datenmodell)
6. [Kantengenerierungs-Engine v3.1](#6-kantengenerierungs-engine-v31)
7. [Werkzeugreferenz](#7-werkzeugreferenz)
8. [Automatisierte Überwachung](#8-automatisierte-überwachung)
9. [Fehlerbibliothek & Problembehandlung](#9-fehlerbibliothek--problembehandlung)
10. [Plattformspezifische Hinweise](#10-plattformspezifische-hinweise)
11. [Bewährte Praktiken](#11-bewährte-praktiken)
12. [FAQ](#12-faq)
13. [Glossar](#13-glossar)

---

## 1. Überblick & Designphilosophie

### 1.1 Was ist das?

Die Literature Learning Suite (LLS) ist ein vollständiges System zur Entdeckung wissenschaftlicher
Literatur, Tiefenanalyse, Wissensgraph-Konstruktion und automatisierten Überwachung. Sie ist kein
„weiteres Literaturverwaltungswerkzeug" — sie ist ein **Literatur-Kognitionsbetriebssystem**.

Herkömmliche Literaturverwaltungsprogramme (Zotero, EndNote, Mendeley) lösen das Problem des
„Wo speichern?". LLS löst die Probleme des **„Wie lesen?"** und **„Wie denken?"**. Durch ein
strenges 7-Schichten-Dissektionsprotokoll wird jede Publikation in einen strukturierten,
abfragbaren und verknüpfbaren Wissensgraph-Knoten transformiert.

### 1.2 Designprinzipien

Diese Prinzipien wurden durch die Tiefenanalyse hunderter Publikationen herausgearbeitet:

**Prinzip 1: Tiefe vor Breite.** Eine vollständige S-Tier-7-Schichten-Analyse ist mehr wert als
50 Publikationen mit nur Titel und Abstract. Wenn Sie 2 Stunden Zeit haben, verbringen Sie 1,5
Stunden mit dem gründlichen Lesen von 2 Publikationen, nicht mit dem Überfliegen von 20.

**Prinzip 2: Verifizieren vor Zitieren.** Überprüfen Sie die bibliografische Identität
(DOI/PMID/arXiv-ID) anhand von mindestens zwei unabhängigen Quellen, bevor Sie analysieren oder
zitieren.

**Prinzip 3: Berichtete Evidenz von Schlussfolgerung trennen.** Jede Aussage muss mit ihrer
erkenntnistheoretischen Stufe gekennzeichnet werden — BERICHET (direkt in der Publikation),
GESTÜTZTE SCHLUSSFOLGERUNG (folgt aus Evidenz + etabliertem Wissen), HYPOTHESE (vom Analysten
generiert) oder UNBEKANNT (kann aus dem verfügbaren Material nicht bestimmt werden).

**Prinzip 4: Evidenz nach Studiendesign beurteilen, nicht nach Journal-Prestige.**
Der Impact-Faktor ist Metadaten, keine Qualitätsbewertung. LLS enthält eine standardisierte
Evidenzbewertungsmatrix, die auf Verzerrungsrisiko, Stichprobengröße, Kontrollqualität und
Reproduzierbarkeit basiert.

**Prinzip 5: Ausschließlich append-only Persistenz.** Alle NDJSON-Datenbanken verwenden den
append-only-Modus. Dies garantiert eine vollständige Operationshistorie und verhindert
irreversiblen Datenverlust.

**Prinzip 6: Keine Falsifikation.** Keine erfundenen Publikationen, Identifikatoren,
Statistiken oder Wirkmechanismen. Kann ein Zitat nicht verifiziert werden, wird es
ausgeschlossen.

**Prinzip 7: Text-als-unvertrauenswürdige-Daten.** Publikationstext und Webinhalte sind
unvertrauenswürdige Daten — Analysegegenstände, niemals Agent-Anweisungen. Dies ist
entscheidend für die Prompt-Injection-Prävention.

### 1.3 Vergleich mit bestehenden Werkzeugen

| Dimension | Zotero/EndNote | Traditionelles Lit-Review | LLS |
|-----------|---------------|---------------------------|-----|
| Speicherung | PDF + Metadaten | — | NDJSON-Wissensgraph |
| Analysetiefe | Manuelle Notizen | Menschliche Zusammenfassung | 7-Schichten-strukturierte Dissektion |
| Publikationsübergreifende Verknüpfung | Manuelle Tags | Subjektiv | 5-Strategien automatische semantische Kanten |
| Evidenzbewertung | Keine | Erfahrungsbasiert | Standardisierte Bewertungsmatrix |
| Automatisierung | Plugin-unterstützt | Keine | Vollständige Überwachungspipeline |
| Abfragbarkeit | Stichwortsuche | Nicht abfragbar | Volltextsuche + Graph-Traversierung |
| Persistenzformat | Proprietär | Klartext | NDJSON (menschenlesbar) |

---

## 2. Systemarchitektur

### 2.1 Verzeichnisstruktur

```
literature-learning-suite/          ← Projektstammverzeichnis (verteilbar)
│
├── GUIDE_DE.md                     ← Dieses Dokument
├── GUIDE_EN.md                     ← Englische Anleitung
├── GUIDE_ZH.md                     ← Chinesische Anleitung (primär)
├── GUIDE_JA.md                     ← Japanische Anleitung
├── GUIDE_KO.md                     ← Koreanische Anleitung
├── SKILL.md                        ← Agent-Betriebsprotokoll
├── README.md                       ← Projektübersicht
├── LICENSE                         ← CC BY-NC-SA 4.0
│
├── scripts/                        ← Kern-Werkzeugkette (23+ Python + Node.js + R)
│   ├── init_workspace.py           ← Arbeitsbereich-Initialisierung
│   ├── literature_search.py        ← Multi-Quellen-Recherche
│   ├── search_arxiv.py             ← arXiv-spezifische Suche
│   ├── download_biorxiv_api.py     ← bioRxiv/medRxiv API-Bulk-Download
│   ├── verify_citation.py          ← Identitäts-Kreuzvalidierung
│   ├── normalize_records.py        ← Standardisierung & Deduplizierung
│   ├── validate_records.py         ← JSON-Schema-Validierung
│   ├── fulltext_fetch.py           ← Vereinter Volltext-Downloader
│   ├── extract_pymupdf.py          ← PDF Text-/Tabellenextraktion
│   ├── extract_marker.py           ← PDF OCR-Extraktion
│   ├── extract_biorxiv_cdp.mjs     ← Chrome CDP-Extraktion (Node.js)
│   ├── biorxiv_chrome_cdp_launcher.bat ← CDP-Launcher (Windows)
│   ├── kg.py                       ← KG-CLI: hinzufügen/statistik/suchen/audit
│   ├── kg_core.py                  ← KG-Kernbibliothek: CRUD + automatische IF-Annotation
│   ├── ll_common.py                ← Gemeinsame Werkzeuge
│   ├── workspace_paths.py          ← Laufzeit-Pfadauflösung
│   ├── gen_edges.py                ← Semantische Kantengenerierung v3.1
│   ├── gen_digest.py               ← Tagesbericht-Generator
│   ├── build_network.py            ← Interaktiver kräftebasierter Netzwerkgraph
│   ├── selfcheck_knowledge_graph.py ← 10-Dimensionen-Qualitätsaudit
│   ├── export_citations.py         ← BibTeX/CSL-Export
│   ├── export_bioc_genes.R         ← Gen-/Pathway-Wörterbuch-Regenerator
│   ├── journal_metrics.py          ← Journal-Metriken-Import
│   ├── monitor.py                  ← Fortsetzbarer Batch-Monitor
│   ├── check_assets.py             ← Asset-Integritätsprüfung
│   └── requirements.txt            ← Python-Abhängigkeiten
│
├── assets/
│   ├── data/                       ← Validierte Referenzdaten (gebündelt)
│   │   ├── bioc_genes.json         ← 90.125 Human- + Maus-Gensymbole
│   │   ├── kegg_pathways.json      ← 25.939 KEGG-Pathways + GO-BP-Terme
│   │   ├── journal_metrics_2024.json ← 21.800 Journal-IF/Quartile
│   │   ├── data-manifest.json      ← SHA-256-Prüfsummen + Herkunftsnachweis
│   │   ├── evidence-rubric.json    ← Evidenzbewertungsmatrix
│   │   ├── relation-ontology.json   ← Kanten-Ontologie
│   │   ├── study-designs.json      ← Studiendesign-Klassifikation
│   │   ├── search-query-packs.json ← Vorkonfigurierte Suchvorlagen
│   │   └── arxiv-categories.json   ← arXiv-Kategoriensystem
│   ├── schemas/                    ← JSON-Schemata (7 Datensatztypen)
│   └── templates/                  ← Datensatz- und Konfigurationsvorlagen
│
├── references/                     ← Methodik & Protokolle (25+ Dokumente)
│   ├── data-model.md               ← Datenmodell-Spezifikation
│   ├── deep-analysis-protocol.md   ← 7-Schichten-Analyse-Detailmethodik
│   ├── s-tier-audit.md             ← Leerhüllen-S-Erkennungsregeln
│   ├── s-tier-examples.md          ← Analysebeispiele nach Publikationstyp
│   ├── s-tier-upgrade-workflow.md  ← Massen-Upgrade-Workflow
│   ├── llm-deep-reasoning-examples.md ← LLM-Schlussfolgerungsmuster
│   ├── gen-edges-v3.md             ← Kantenalgorithmus-Detail
│   ├── edge-generation.md          ← Kantentyp-Referenz
│   ├── bioconductor-entity-matching.md ← Gen-/Pathway-Abgleichslogik
│   ├── full-text-access.md         ← Rechtskonforme Volltextbeschaffung
│   ├── preprint-fulltext.md        ← bioRxiv/medRxiv-Extraktion
│   ├── self-review-checklist.md    ← Checkliste nach Code-Änderungen
│   ├── connectivity.md             ← Netzwerk-/Proxy-Konfiguration
│   ├── cron-troubleshooting.md     ← Unbeaufsichtigte Ausführungsdiagnose
│   ├── automation.md               ← Überwachungsautomatisierung
│   ├── hermes-monitoring-template.md ← Agent-Cron-Vorlage
│   ├── mcp-integration.md          ← MCP-Integrationsleitfaden
│   ├── mcp-and-tool-routing.md     ← Werkzeug-Fallback-Strategie
│   ├── journal-metrics-2024.md     ← JCR-Metriken-Nutzungshinweise
│   ├── pdf-and-ocr.md              ← PDF-Extraktionsmethoden-Vergleich
│   └── bioinfo-tools.md            ← Bioinformatik-Hilfswerkzeuge
│
└── tests/                          ← Unit-Tests + synthetische Testdaten
```

### 2.2 Datenflusspipeline

```
                    Forschungsfrage
                          │
                          ▼
                 ┌─────────────────┐
                 │ 1. Suchstrategie │  ← PICO/PECO-Framework
                 │    definieren    │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 2. Multi-Quellen │  ← PubMed/arXiv/bioRxiv/Crossref
                 │    Recherche     │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 3. Identität     │  ← Zweiquellen-Kreuzvalidierung DOI/PMID
                 │    verifizieren  │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 4. Volltext      │  ← PMC XML / arXiv HTML / CDP-Extraktion
                 │    beschaffen    │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 5. S-Tier-7-    │  ← LLM-Tiefenanalyse (T1–T7)
                 │    Schichten     │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 6. Evidenz       │  ← Verzerrungsrisiko/Stichprobengröße/
                 │    bewerten      │     Reproduzierbarkeit
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ 7. Wissensgraph  │  ← NDJSON append
                 │  ┌───────────┐  │
                 │  │ papers.db │  │
                 │  │concepts.db│  │
                 │  │ edges.db  │  │
                 │  └───────────┘  │
                 └────────┬────────┘
                          │
            ┌─────────────┼─────────────┐
            ▼             ▼             ▼
      ┌──────────┐ ┌──────────┐ ┌──────────┐
      │  Kanten  │ │ Tages-   │ │ Netzwerk │
      │gen_edges │ │ bericht  │ │build_net │
      └──────────┘ └──────────┘ └──────────┘
            │             │             │
            ▼             ▼             ▼
      ┌──────────────────────────────────┐
      │     9. Qualitätsselbstkontrolle   │
      │        + Automat. Überwachung     │
      └──────────────────────────────────┘
```

### 2.3 Laufzeit-Arbeitsbereich

Der durch `init_workspace.py` erstellte Arbeitsbereich ist ein eigenständiges Verzeichnis:

```
my-workspace/                      ← gitignored, zur Laufzeit generiert
├── workspace.json                 ← Arbeitsbereich-Metadaten
├── papers.db                      ← NDJSON, eine Publikation pro Zeile
├── concepts.db                    ← NDJSON, ein Konzeptknoten pro Zeile
├── edges.db                       ← NDJSON, eine semantische Kante pro Zeile
├── queries.db                     ← NDJSON, Suchprotokoll
├── journal_metrics.db             ← NDJSON (optional, benutzerdefinierte IF-Daten)
├── data/                          ← Kopie der eingepflanzten Gen-/Pathway-Wörterbücher
├── fulltext/                      ← Heruntergeladene Volltextdokumente
├── fulltext_cache/                ← Volltext-Cache
├── reports/                       ← Generierte Berichte
├── daily_digest/                  ← Tägliche Zusammenfassungen
├── config/                        ← Überwachungskonfiguration
├── exports/                       ← BibTeX-Exporte
├── imports/                       ← Zu importierende Daten
├── cache/                         ← Allgemeiner Cache
├── biorxiv_api/                   ← bioRxiv-API-Antwort-Cache
└── logs/                          ← Laufzeitprotokolle
```

---

## 3. Schnellstart

### 3.1 Voraussetzungen

- Python 3.10+
- pip
- Node.js v24+ (nur für bioRxiv CDP-Extraktion erforderlich)
- R + Bioconductor (nur für Gen-Wörterbuch-Regeneration erforderlich)

### 3.2 Installation

```bash
# Repository klonen
git clone <repo-url>
cd literature-learning-suite

# Python-Abhängigkeiten installieren
pip install -r scripts/requirements.txt
```

### 3.3 Arbeitsbereich initialisieren

```bash
python scripts/init_workspace.py --root ./my-workspace
```

Dieser Befehl:
1. Erstellt die Arbeitsbereich-Verzeichnisstruktur
2. Kopiert das Gen-Wörterbuch (90.125 Gensymbole) nach `my-workspace/data/`
3. Kopiert Pathway-/GO-Terme (25.939 Labels) nach `my-workspace/data/`
4. Erstellt leere NDJSON-Datenbankdateien (papers/concepts/edges/queries/journal_metrics.db)
5. Installiert die Standard-Überwachungskonfiguration aus Vorlagen

### 3.4 Umgebungsvariable setzen

```bash
# Linux/macOS
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"

# Windows PowerShell
$env:LITERATURE_KG_ROOT = (Resolve-Path ./my-workspace)
```

Ohne diese Variable verwenden Skripte standardmäßig `./literature-workspace`.

### 3.5 Erste Suche

```bash
# PubMed-Recherche
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10

# arXiv-Recherche
python scripts/literature_search.py arxiv "single cell foundation model" -n 10

# Crossref-Recherche
python scripts/literature_search.py crossref "tumor microenvironment review" -n 10

# bioRxiv-Recherche
python scripts/literature_search.py biorxiv "spatial transcriptomics" -n 10
```

### 3.6 Publikationen einpflegen

```bash
# Suchergebnisse standardisieren
python scripts/normalize_records.py search-results.jsonl -o normalized.jsonl

# Datensatzformat validieren
python scripts/validate_records.py normalized.jsonl

# Zum Wissensgraph hinzufügen (automatische Deduplizierung)
python scripts/kg.py --root ./my-workspace add normalized.jsonl

# Statistik anzeigen
python scripts/kg.py --root ./my-workspace stats
```

---

## 4. Kernmethodik: S-Tier-7-Schichten-Protokoll

Dies ist das Herzstück von LLS. Jede Forschungspublikation durchläuft eine vollständige
7-Schichten-Analyse.

### 4.1 Warum 7 Schichten?

Herkömmliches Lesen endet bei „Abstract lesen → taggen → zwei Sätze schreiben". Diese flache
Verarbeitung führt dazu, dass man:
- die tiefere logische Struktur einer Publikation nicht erfasst,
- implizite Verbindungen zwischen Publikationen nicht entdeckt,
- die tatsächliche Evidenzstärke nicht bewerten kann,
- kein abfragbares strukturiertes Wissen aufbauen kann.

Das 7-Schichten-Protokoll zerlegt den Leseprozess in sieben unabhängige, aber komplementäre
kognitive Dimensionen und zwingt den Analysten (Mensch oder LLM), in jeder Dimension
substanzielle Ergebnisse zu produzieren.

### 4.2 Schicht 1 (T1): Bibliografisches & Studienprofil

**Ziel**: Verifizierbare Identität und Kontext herstellen.

**Inhalt**:
- Vollständiger Titel
- Erstautor, Korrespondenzautor und deren Institutionen
- Journalname (mit automatischer JCR-2024-IF- und Quartil-Annotation via `enrich_paper_if()`)
- Stabile Identifikatoren (PMID / DOI / arXiv-ID)
- Dokumenttyp (Article / Review / Preprint / Clinical Trial usw.)
- Studiendesign (RCT / Kohorte / Fall-Kontrolle / Querschnitt usw.)
- Probensystem/Modellsystem, Stichprobengröße
- Volltextstatus (fulltext / abstract_only / metadata_only / unavailable)
- Volltextquelle, Zugriffsdatum, Version

**Automatisierungsgrad**: Diese Schicht kann werkzeuggestützt automatisch extrahiert werden.

### 4.3 Schicht 2 (T2): Kernfragestellung

**Ziel**: Die Publikation in **eine falsifizierbare Kernfrage** destillieren und in **≥5
überprüfbare Teilfragen** zerlegen.

**Anforderung**: Die Kernfrage muss **mechanistisch** sein, nicht deskriptiv. Nicht „Was haben
sie entdeckt?", sondern „Wie verursacht X durch Y das Phänomen Z?".

**Beispiel (gut)**:
> Kernfrage: Wie induziert GDF15 über Xenobiotika-Rezeptor-Signalwege eine NK-Zell-Dysfunktion
> in der Tumor-Mikroumgebung?
>
> Teilfragen:
> 1. In welchen Tumortypen wird GDF15 überexprimiert? Zusammenhang mit der Prognose?
> 2. Welche Rezeptoren vermitteln das GDF15-Signal? Sind nicht-kanonische Rezeptoren beteiligt?
> 3. Wie verändert die GDF15-Signalkaskade die zytotoxische Funktion von NK-Zellen?
> 4. Ist diese Dysfunktion reversibel? Wo liegen pharmakologische Angriffspunkte?
> 5. Stellt die Blockade des GDF15-Signals die Anti-Tumor-Immunität in vivo wieder her?

**Beispiel (schlecht — Leerhülle S)**:
> Kernfrage: Diese Arbeit untersuchte GDF15.
> Teilfragen: (leer)

### 4.4 Schicht 3 (T3): Behauptung-Evidenz-Synthese-Ketten (≥5)

**Ziel**: Die Kernbefunde in ≥5 unabhängige C-E-S-Ketten mit konkreter Evidenz transformieren.

**Struktur jeder Kette**:

| Feld | Anforderung | Schwellenwert |
|------|------------|---------------|
| `claim` | Falsifizierbare Schlussfolgerung in eigenen Worten | — |
| `evidence` | Konkrete Daten: N, Effektstärke, p-Wert, Assay | **>20 Zeichen** (Leerhüllen-Trigger) |
| `synthesis` | Wie stützt/schwächt die Evidenz die Behauptung? Domänenübergreifende Verknüpfungen? | — |
| `strength` | 1–5★, basierend auf Evidenzqualität | — |
| `uncertain` | Alternativerklärungen, fehlende Validierung, Störfaktoren | — |

**Beispiel (gut)**:
> Claim: Die PD-L1-1–49%-Subgruppe profitiert von neoadjuvanter Chemoimmuntherapie.
> Evidence: 2847 NSCLC-Patienten, SHAP-Interaktionsanalyse, EFS HR=0,52 (95% CI 0,34–0,78),
> p=0,002, unabhängige radiologische Begutachtung.
> Synthesis: Schließt die Evidenzlücke für die PD-L1-„Grauzonen"-Population, die in der
> KEYNOTE-671-Subgruppenanalyse ausgeschlossen war, und deutet darauf hin, dass das
> CheckMate-816-Regime auf eine breitere Population anwendbar sein könnte.
> Strength: ★★★★ (multizentrische RCT, große Stichprobe, unabhängige Begutachtung)
> Uncertain: Inter-Assay-Konkordanz der PD-L1-Testplattformen (22C3 vs. 28-8); hoher Anteil
> asiatischer Patienten könnte Übertragbarkeit einschränken.

**Beispiel (schlecht — Leerhülle S)**:
> Claim: Ein neuer Biomarker wurde entdeckt.
> Evidence: Die Autoren bewiesen ihre Hypothese.
> (evidence-Feld ≤20 Zeichen → Leerhüllen-Erkennung schlägt an)

### 4.5 Schicht 4 (T4): Mechanismus-Kaskade

**Ziel**: Die vollständige Kausalkette vom Trigger zum Phänotyp abbilden, präzise bis zur
Modifikationsstelle.

**Struktur**:
```
Trigger
  │
  ▼
[Upstream-Rezeptor] → [Second Messenger/Kinase] → [Transkriptionsfaktor] → [Zielgen] → [Phänotyp]
```

**Muss enthalten**:
- **≥3 Kausalschritte** mit Direktionalität
- **Schlüssel-Modifikationsstellen**: präzise bis zum Aminosäurerest (z.B. „NF-κB p65 Ser536
  Phosphorylierung", nicht „NF-κB-Aktivierung")
- **Downstream-Effekte**: Zellverhalten, Stoffwechsel, Interaktionsänderungen
- **Feedback-Schleife**: ≥1 positive oder negative Rückkopplung

**Evidenzstatus-Annotation** (für jeden Schritt erforderlich):
- `demonstrated` (demonstriert): in dieser Arbeit direkt validiert
- `supported` (gestützt): durch indirekte Evidenz gestützt
- `background` (Hintergrund): allgemein anerkanntes Wissen
- `hypothesis` (Hypothese): Vermutung des Analysten

**Beispiel**:
```
Trigger: Tumor-abgeleitetes GDF15
  │
  ▼
GFRAL-Rezeptorbindung → [demonstriert: Co-IP, Abb. 2A]
  │
  ▼
JAK2-STAT3-Signalwegaktivierung → [gestützt: Phospho-Antikörper, Abb. 3B]
  │
  ▼
STAT3 pTyr705-Phosphorylierung + Kerntranslokation → [demonstriert: subzelluläre Fraktionierung+WB, Abb. 3C]
  │
  ▼
SOCS3-Transkriptionsaktivierung (negative Rückkopplung) → [demonstriert: qPCR+Promoter-Luciferase]
  │
  ▼
NK-Zell-Zytotoxizität vermindert (CD107a↓, IFN-γ↓) → [demonstriert: Durchflusszytometrie, Abb. 4]
```

**Verboten**:
- ❌ „Aktiviert Signalweg" (zu vage)
- ❌ Modifikationsstellen frei erfinden
- ❌ Evidenzstatus nicht annotieren

### 4.6 Schicht 5 (T5): Verborgene Organisationsachsen (≥3)

**Ziel**: Muster aufdecken, die die Publikation **nicht explizit** benennt — implizite Annahmen,
räumlich-zeitliche Organisation, Selektionseffekte oder die tiefere experimentelle Logik.

**Jede Achse enthält**:
- `observation`: eine verifizierbare Tatsache aus der Publikation
- `interpretation`: das vom Analysten entdeckte tiefere Muster

**Diese Schicht prüft Ihre Synthesefähigkeit**, nicht Ihre Fähigkeit, den Diskussionsteil zu
kopieren.

**Beispiel**:
> Observation 1: In allen Dimensionsreduktionsanalysen clustern Tumorrand-Proben stets getrennt
> von Tumorzentrum-Proben.
> Interpretation 1: Die Studie definiert implizit den „Rand" und nicht das „Zentrum" als das
> krankheitsbestimmende Kompartiment. Dies erklärt, warum Genexpressions-Signaturen des Zentrums
> eine geringere prognostische Aussagekraft haben — die molekulare Heterogenität des Zentrums
> wird durch das Mikroumgebungssignal des Randes überlagert.
>
> Observation 2: In der Einzelzellanalyse treten Veränderungen in Immun-Subpopulationen vor
> Veränderungen in Tumor-Subpopulationen auf.
> Interpretation 2: Die Zeitreihendaten deuten auf ein „Immun-zuerst"-Krankheitsprogressionsmodell
> hin — die Mikroumgebungs-Remodellierung ist Treiber, nicht Folge der Tumorevolution. Trifft
> diese Interpretation zu, sollten frühe Interventionsziele auf Immunzellen, nicht auf
> Tumorzellen liegen.
>
> Observation 3: Unter den differentiell exprimierten Genen zwischen Respondern und
> Non-Respondern ist die Anreicherung von Stoffwechselwegen (45%) deutlich höher als die von
> Immunwegen (12%).
> Interpretation 3: Obwohl Titel und Diskussion der Publikation auf Immunmechanismen fokussieren,
> zeigt die intrinsische Datenstruktur, dass metabolische Reprogrammierung der fundamentalere
> Bestimmungsfaktor sein könnte. Es besteht eine systematische Verschiebung zwischen dem
> Narrativ der Publikation (Immun→Wirksamkeit) und ihrer Datengewichtung (Metabolismus→Wirksamkeit).

**Verboten**:
- ❌ T5 = Diskussionsteil der Publikation paraphrasieren
- ❌ T5 = „Die Autoren entdeckten X" (das ist T3, nicht T5)

### 4.7 Schicht 6 (T6): Konzeptueller Beitrag

**Ziel**: Den substanziellen wissenschaftlichen Beitrag identifizieren, nicht abstrakt
„untersuchte die Rolle von X in Y".

**Vier Elemente**:

1. **Neue Konzepte** (≥1): Mit operationaler Definition, unabhängig zitierbar.
2. **Widerlegte/revidierte Ansichten** (≥1): Welche etablierten Überzeugungen werden
   herausgefordert oder eingeschränkt?
3. **Methodologische Durchbrüche** (≥1): Welche neue technische Fähigkeit, die andere nutzen
   können, wird beigetragen?
4. **Randbedingungen**: Unter welchen Bedingungen gilt dieser Beitrag **nicht**?

**Beispiel**:
> Neues Konzept:
> - „Immunmetabolischer Checkpoint" (Immunometabolic Checkpoint): Definiert als
>   metabolitvermittelte Immunsuppressionsachse über Immunrezeptoren, abzugrenzen von
>   klassischen Protein-Ligand-Rezeptor-Immuncheckpoints.
>   - Operationale Definition: Erfordert gleichzeitig (a) niedermolekularen Metaboliten als
>     Liganden, (b) Signaltransduktion über Immunrezeptoren, (c) reversible Immunzell-
>     Funktionssuppression.
>
> Widerlegte Ansicht:
> - Revidiert die monofunktionale Annahme „GDF15 ist nur ein Appetitzügler" und belegt eine
>   appetitunabhängige Rolle in der Tumorimmunologie.
>
> Methodologischer Durchbruch:
> - Etablierung einer antikörperfreien Metabolit-Rezeptor-Bindungsdetektionsmethode
>   (SILAC-Markierung + chemische Quervernetzung + Massenspektrometrie), übertragbar auf
>   andere orphan-Metabolitrezeptoren.
>
> Randbedingungen:
> - Mechanismus gilt nur in Tumoren mit hoher GDF15-Expression (≥ 2× Median); in
>   niedrig-exprimierenden Tumoren ist diese Achse inaktiv.
> - Nicht in immundefizienten Modellen validiert; synergistische Wirkung der adaptiven
>   Immunität nicht ausgeschlossen.

### 4.8 Schicht 7 (T7): Publikationsübergreifende Beziehungen (≥5)

**Ziel**: Die Publikation mittels substanzieller biologischer Beziehungen mit anderen
verifizierten Datensätzen verknüpfen.

**Anforderungen pro Beziehung**:
- `ref_id`: Stabiler Identifikator der Zielpublikation (PMID:xxxxx oder DOI:10.xxxx/xxxxx)
- `relation`: Beziehungstyp (siehe gültige Typen unten)
- `description`: **60–150 Wörter** Erklärung, WARUM diese Beziehung besteht

**Gültige Beziehungstypen**:
- `supports` (stützt): Evidenz dieser Arbeit stützt die Schlussfolgerung der Zielpublikation
- `contradicts` (widerspricht): Evidenz dieser Arbeit widerspricht der Zielpublikation
- `extends` (erweitert): Diese Arbeit erweitert die Zielpublikation sinnvoll
- `replicates` (repliziert): Unabhängige Reproduktion der Kernbefunde der Zielpublikation
- `methodological_complement` (methodische Ergänzung): Gleiche Frage, andere Technik
- `shared_mechanism` (gemeinsamer Mechanismus): Geteilter molekularer Mechanismus
- `upstream_of` / `downstream_of` (oberhalb/unterhalb in der Kausalkette)
- `clinical_translation` (klinische Translation): Grundlagenforschung zu Klinik
- `shares_disease_model` (gemeinsames Krankheitsmodell)

**Ausdrücklich verbotene nicht-biologische Beziehungen** (werden von `gen_edges.py` Strategie 1
herausgefiltert):
- `same_journal` (gleiches Journal)
- `same_issue` (gleiche Ausgabe)
- `same_author` (gleicher Autor)
- `same_year` (gleiches Jahr)

### 4.9 Leerhüllen-S-Erkennung

Ein als S-Tier markierter Datensatz ist eine Leerhülle, wenn **irgendeine** der folgenden
Bedingungen zutrifft:

1. `tier2_subquestions` ist leer oder fehlt
2. `tier3_ces_chains` enthält weniger als 5 Ketten
3. Irgendein `evidence`-Feld ≤ 20 Zeichen
4. `tier4_mechanism_cascade` hat weniger als 3 Kaskadenschritte
5. `tier5_hidden_axis` enthält weniger als 3 Achsen
6. `tier7_cross_refs` enthält weniger als 5 Einträge
7. Irgendeine T7-Beziehung gehört zur Verbotsliste

**Post-Write-Verifikation**:
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
              f'T3={len(chains)} empty_evidence={empty_ev} '
              f'T4_steps={len(p.get(\"tier4_mechanism_cascade\",{}).get(\"cascade\",[]))} '
              f'T5={len(p.get(\"tier5_hidden_axis\",[]))} '
              f'T7={len(p.get(\"tier7_cross_refs\",[]))}')
"
```

### 4.10 Nicht-Forschungspublikationen

Bei Publikationen der folgenden Typen wird keine S-Tier-Analyse durchgeführt:
- **Nachrichten/Editorials** (news/editorial)
- **Errata/Zurückziehungen** (erratum/retraction)
- **Research Briefs/Communications** ohne substanzielle Methoden-/Ergebnisinhalte

Für diese Publikationen: `analysis_tier: "NR"` (Non-Research) setzen, den tatsächlichen Typ in
`core_findings` eintragen und **keine Analyse-Inhalte erfinden**.

---

## 5. Wissensgraph-Datenmodell

### 5.1 Designprinzipien

LLS verwendet NDJSON (Newline-Delimited JSON) anstelle von SQLite. Warum?

1. **Menschenlesbar und -editierbar**: Jede Zeile ist ein vollständiges JSON-Objekt
2. **Versionskontrollfreundlich**: `git diff` zeigt Änderungen zeilenweise
3. **Append-only-Datenintegrität**: Unveränderlicher Verlauf, kein irreversibler Datenverlust
4. **CLI-abfragbar**: Standardwerkzeuge wie `head`, `grep`, `jq` direkt nutzbar
5. **Keine externen Abhängigkeiten**: Kein Datenbanktreiber erforderlich

### 5.2 papers.db

Hauptdatenbank der Publikationen, ein vollständiger Datensatz pro Zeile.

**Pflichtfelder**: `id` (stabiler Identifikator), `title`, `source` (Herkunft), `retrieved_at`

**S-Tier-Analysefelder** (v4.0-Standard):
```json
{
  "id": "PMID:42251595",
  "title": "Vollständiger Titel",
  "authors": ["Erstautor", "..."],
  "journal": "Vollständiger Journalname",
  "impact_factor": 63.1,
  "journal_quartile": "Q1",
  "doi": "10.xxxx/xxxxx",
  "pmid": "42251595",
  "source": "pubmed",
  "retrieved_at": "2026-06-09T09:00:00",
  "fulltext_status": "fulltext",
  "analysis_tier": "S",
  "analysis_method": "LLM_deep_reasoning_S_tier_v4.0",
  "tier2_core_question": "...",
  "tier2_subquestions": ["Q1", "Q2", "Q3", "Q4", "Q5"],
  "tier3_ces_chains": [
    {"chain_id":1, "claim":"...", "evidence":"...", "synthesis":"...", "strength":3, "uncertain":"..."}
  ],
  "tier4_mechanism_cascade": {
    "trigger":"...",
    "cascade":["Schritt1","Schritt2","Schritt3"],
    "key_modifications":[...],
    "downstream_effects":"...",
    "feedback":[...],
    "evidential_status":{...}
  },
  "tier5_hidden_axis": [{"observation":"...", "interpretation":"..."}],
  "tier6_concept_innovation": {
    "new_concepts":[...],
    "overturned_views":[...],
    "methodological_breakthroughs":[...],
    "boundary_conditions":"..."
  },
  "tier7_cross_refs": [{"ref_id":"PMID:xxxxx", "relation":"extends", "description":"..."}],
  "entities": {
    "genes":["GDF15"],
    "pathways":["JAK-STAT"],
    "cell_types":["NK-Zellen"],
    "diseases":["Krebs"]
  }
}
```

### 5.3 concepts.db

Konzept-/Entitätsdatenbank.

```json
{
  "id": "CONCEPT:immunmetabolischer_checkpoint",
  "name": "Immunmetabolischer Checkpoint",
  "type": "mechanism",
  "definition": "Metabolit-vermittelte Immunsuppressionsachse über Immunrezeptoren",
  "source_papers": ["PMID:42251595", "PMID:39988000"],
  "created_at": "2026-06-09T09:00:00"
}
```

Gültige `type`-Werte: `mechanism`, `disease`, `method`, `cell_type`, `pathway`, `drug`, `gene`,
`phenomenon`, `hypothesis`

### 5.4 edges.db

Datenbank der semantischen Kanten.

```json
{
  "source": "PMID:42251595",
  "target": "PMID:39988000",
  "relation": "extends",
  "description": "PMID:42251595 erweitert die in PMID:39988000 identifizierte GDF15-GFRAL-Bindung um die Entdeckung der JAK2-STAT3-SOCS3-Negativrückkopplungsschleife und erweitert damit das GDF15-Signal von der Ligand-Rezeptor-Erkennung auf eine vollständige Signaltransduktionskaskade. Beide Arbeiten bilden gemeinsam den molekularen Rahmen der GDF15-Immunsuppressionsachse.",
  "provenance": "analyst",
  "created_at": "2026-06-09T10:00:00"
}
```

### 5.5 queries.db

Suchprotokoll für Audits und Reproduzierbarkeit.

```json
{
  "source": "pubmed",
  "query": "spatial transcriptomics cancer",
  "executed_at": "2026-06-09T09:00:00",
  "result_count": 42,
  "parameters": {"retmax": 50, "sort": "date"}
}
```

---

## 6. Kantengenerierungs-Engine v3.1

### 6.1 Überblick

`gen_edges.py` erstellt biologisch bedeutsame semantische Verbindungen zwischen Publikationen
mittels 5 komplementärer Strategien, alle implementiert mit invertierten Indizes (O(M)-Komplexität
statt O(N²)).

### 6.2 Fünf Strategien

| # | Strategie | Relation | Logik |
|---|-----------|----------|-------|
| 1 | Explizite Referenzen | extends, accompanied_by | tier7_cross_refs aus der Tiefenanalyse (nicht-bio gefiltert) |
| 2 | Gemeinsame Moleküle | shares_molecules (≥2) | 90.125 Gene, zweistufiger Abgleich |
| 2.5 | Textüberlappung | shares_topic (≥4) | Bag-of-Words aus Kernbefunden → invertierter Index |
| 3 | Erkrankung × Methode | shares_disease_method | Kreuzprodukt Krankheitslabel × Methodenlabel |
| 4 | Verborgene Achse | shares_paradigm | T5-Tiefenmuster-Schlüsselwortresonanz |
| 5 | Konzeptknoten | defines_concept | Publikation→Konzept-Kanten aus concepts.db |

#### Strategie 1: Explizite Referenzen
- **Datenquelle**: `tier7_cross_refs`-Feld (Ergebnis der LLM-Tiefenanalyse)
- **Filter**: Ausschluss nicht-biologischer Beziehungen
- **Qualität**: Höchste (menschliches/LLM-Expertenurteil), aber abhängig von Analysetiefe

#### Strategie 2: Gemeinsame Moleküle (≥2 gemeinsam)
- **Datenquelle**: 90.125 Human- und Maus-Gensymbole (Bioconductor-Export)
- **Abgleichslogik**: Zweistufig — (1) Regex-Extraktion von Kandidaten aus Text, (2) O(1)-Lookup
  im Genset zur Bestätigung
- **Vorteil**: Präzise (eliminiert False-Positives von ML/TNF-Abkürzungen)
- **Limitierung**: Maximal 15 Publikationen pro Gen (verhindert Kantenexplosion bei
  hochfrequenten Genen wie TNF)

#### Strategie 2.5: Textüberlappung (≥4 gemeinsame Schlüsselwörter)
- **Datenquelle**: Kernbefunde (`core_findings` oder T3-Claims)
- **Abgleichslogik**: (1) Bag-of-Words (ohne Stoppwörter), (2) Wort→Publikation invertierter
  Index, (3) Zählung gemeinsamer Wörter pro Publikationspaar
- **Anteil**: Hauptstrategie, liefert ca. 70% aller Kanten

#### Strategie 3: Gleiche Erkrankung, gleiche Methode
- **Datenquelle**: `diseases`- und `methods`/`technologies`-Felder
- **Abgleichslogik**: Kreuzprodukt Krankheitslabel × Methodenlabel
- **Charakteristik**: Größere Granularität, niedrigere Abdeckung, aber hohe Präzision

#### Strategie 4: Verborgene Achsen
- **Datenquelle**: `tier5_hidden_axis`-Feld
- **Abgleichslogik**: Extraktion von Tiefenmuster-Schlüsselwörtern (Paradigma, Bias, Survivor,
  Selektion) und Resonanzabgleich unter den Top-200-Publikationen
- **Charakteristik**: Tiefste Kantenart, erfasst implizite Paradigmen-Gemeinsamkeiten

#### Strategie 5: Konzeptknoten
- **Datenquelle**: `concepts.db`
- **Abgleichslogik**: Für jedes Konzept `defines_concept`-Kanten zu seinen `source_papers`
  generieren
- **Charakteristik**: Verbindet Konzeptknoten mit Publikationsknoten

### 6.3 Verwendung

```bash
# KRITISCH: Bytecode-Cache vor jedem Lauf leeren
rm -rf scripts/__pycache__

# Kantengenerator ausführen (-B verhindert Schreiben neuer .pyc)
python -B scripts/gen_edges.py

# Interaktiven Netzwerkgraph aktualisieren
python scripts/build_network.py
```

### 6.4 Ausgabeformat

```json
{
  "source": "PMID:42251595",
  "target": "PMID:39988000",
  "relation": "shares_molecules",
  "description": "Gemeinsame Moleküle: genes:GDF15, genes:TNF, pathways:JAK-STAT signaling",
  "metadata": {"shared_entities": ["genes:GDF15", "genes:TNF", "pathways:JAK-STAT signaling"]},
  "provenance": "deterministic_molecule_index",
  "created_at": "2026-06-09T10:00:00"
}
```

---

## 7. Werkzeugreferenz

### 7.1 Recherche & Entdeckung

#### `literature_search.py`
Multi-Quellen-Recherche in wissenschaftlichen Datenbanken.

```
Verwendung: literature_search.py <quelle> <suchbegriff> [optionen]

Quellen (source):
  pubmed    - PubMed E-utilities API
  arxiv     - arXiv API
  crossref  - Crossref API
  biorxiv   - bioRxiv API

Optionen:
  -n N      - Maximale Trefferanzahl (Standard: 20)
  -y YEAR   - Nach Jahr filtern
  -o FILE   - Ausgabedatei (Standard: stdout)
```

#### `search_arxiv.py`
arXiv-spezifische Suche mit Filterung nach Autor, Kategorie und Datumsbereich.

```
Verwendung: search_arxiv.py [--author NAME] [--category CAT] [--max N]
```

#### `download_biorxiv_api.py`
bioRxiv/medRxiv Bulk-Download von Abstracts nach Datum.

```
Verwendung: download_biorxiv_api.py --date YYYY-MM-DD [--source biorxiv|medrxiv]
```

### 7.2 Verifikation & Standardisierung

#### `verify_citation.py`
Kreuzvalidierung der Publikationsidentität.

```
Verwendung: verify_citation.py --doi 10.xxxx/xxxxx
            verify_citation.py --pmid 12345678
            verify_citation.py --arxiv 2401.01234
```

Validierungsablauf:
1. Existenz des Datensatzes in der Primärquelle bestätigen
2. Übereinstimmung von Titel, Autor, Jahr, Journal prüfen
3. Auf Zurückziehungen/Errata/Bedenkenhinweise prüfen
4. Bei wichtigen Datensätzen: Kreuzvalidierung in einer zweiten Quelle

#### `normalize_records.py`
Standardisierung und Deduplizierung von Suchergebnissen.

Deduplizierungspriorität: PMID > arXiv-ID > normalisierte DOI > normalisierter Titel+Jahr

### 7.3 Volltextbeschaffung

#### `fulltext_fetch.py`
Vereinter Volltext-Downloader mit automatischer Pfadauswahl.

```
Pfadpriorität:
1. PubMed Central OA XML (benötigt PMCID)
2. arXiv HTML (keine Authentifizierung nötig)
3. bioRxiv/medRxiv API-Abstract (immer verfügbar)
4. Vom Benutzer bereitgestelltes lokales PDF
```

#### `extract_biorxiv_cdp.mjs`
Chrome DevTools Protocol Volltext-Extraktor für Cloudflare-geschützte bioRxiv/medRxiv-Seiten.

```
Voraussetzungen:
1. Chrome im Remote-Debugging-Modus starten (Port 9223)
2. Sicherheitsüberprüfung manuell im Browser absolvieren (nicht automatisiert)
3. Dann dieses Skript ausführen, um den gerenderten Text zu extrahieren

Verwendung: node scripts/extract_biorxiv_cdp.mjs --doi 10.1101/XXXX --port 9223
```

**Designprinzip**: Keine automatisierten CAPTCHA-Umgehungen, keine Zugriffskontrollumgehung.
Der CDP-Extraktor ersetzt lediglich den manuellen „Auswählen→Kopieren→Einfügen"-Vorgang.

### 7.4 Wissensgraph-Operationen

#### `kg.py`
Kommandozeilenwerkzeug für den Wissensgraph.

```
Verwendung: kg.py --root <workspace> <befehl> [argumente]

Befehle:
  add <datei>    - Publikationsdatensätze hinzufügen (JSONL/NDJSON), automatische Deduplizierung
  stats          - Statistik anzeigen (Datensatzanzahl, Quellenverteilung, Volltextstatus)
  search <query> - Volltextsuche (AND-Logik, Groß-/Kleinschreibung ignoriert)
  audit          - Integritätsaudit (doppelte IDs, fehlende Pflichtfelder, JSON-Parsefehler)
```

#### `kg_core.py`
Wissensgraph-Kernbibliothek (Python-API).

```python
from kg_core import (
    add_paper,          # Publikation hinzufügen (automatische Dedup+IF-Annotation)
    enrich_paper_if,    # Automatische JCR-2024-IF-Annotation
    get_stats,          # Statistik abrufen
    get_recent_papers,  # Publikationen der letzten N Tage abrufen
    search_papers,      # Publikationen durchsuchen
    lookup_journal_impact_factor,  # Journal-IF nachschlagen
    write_daily_digest, # Tagesbericht schreiben
)

# Beispiel
paper = {'title': '...', 'pmid': '12345', 'journal': 'Nature', 'source': 'pubmed'}
paper = enrich_paper_if(paper)  # Automatisch IF=50.5, Q1 ergänzt
added = add_paper(paper)  # True (neu) oder False (Duplikat)
```

### 7.5 Analyse & Ausgabe

#### `gen_edges.py`
Semantischer Kantengenerator v3.1 (Details siehe Abschnitt 6).

#### `gen_digest.py`
Generator für den täglichen Zusammenfassungsbericht.

```
Generierte Inhalte:
- Gesamtanzahl Publikationen und S/A/B-Tier-Verteilung
- Verteilung der Kantentypen
- S-Tier-Publikationsliste (Top 15 mit Journal und Kernbefunden)

Ausgabe: daily_digest/JJJJ-MM-TT.md
```

#### `build_network.py`
Generator für interaktive HTML-Netzwerkgraphen.

```
Ausgabe: network.html (Doppelklick zum Öffnen)

Visualisierungsregeln:
- Knotenfarbe: PubMed=Grün, arXiv=Rot, bioRxiv=Blau, medRxiv=Hellblau, Konzept=Gelb
- Kantenfarbe: Publikationsübergreifend=Blau, Gemeinsame Gene=Orange,
  Gemeinsame Pathways=Violett, Publikation→Konzept=Gelb
- Knotengröße ∝ Anzahl der Claims
- Weißer Rand = vollständige Tiefenanalyse
- Rein offline, keine externen Abhängigkeiten
```

#### `selfcheck_knowledge_graph.py`
10-Dimensionen-Qualitätsselbstkontrolle.

| Dimension | Prüfinhalt |
|-----------|-----------|
| Dateiinventar | Vollständiger Verzeichnisscan, Dateianzahl/-größe |
| Verbotene Rückstände | chrome_cdp_profile, Cookie-Dateien, temporäre Testdateien |
| CDP-Port | Ob Port 9222/9223 noch geöffnet ist |
| DB-Integrität | NDJSON-Parsefehler, doppelte IDs, fehlende Pflichtfelder |
| S-Tier-Qualität | Leerhüllen-S / Schwaches-S-Erkennung (7 Kriterien pro Publikation) |
| Konzept-Audit | Doppelte IDs, fehlende Namen |
| Kanten-Audit | Illegale nicht-bio-Beziehungen, unbeschriebene Kanten, verwaiste Kanten, Selbstschleifen, Duplikate |
| Volltext-Cache | Namenskonvention, zu kleine Dateien, Cloudflare-Rückstände, Duplikate |
| Kantenstatistik | Kantenanzahl pro Strategie |
| Netzwerkkonsistenz | Gültigkeit der Querverweise |

#### `export_citations.py`
Export in BibTeX/CSL-Format.

```
Verwendung: export_citations.py <papers.db> --format bibtex -o bibliothek.bib
```

---

## 8. Automatisierte Überwachung

### 8.1 Überwachungskonfiguration

Überwachungsjobs werden in `config/monitor-job.json` definiert:

```json
{
  "name": "Tägliche Top-Journal-Literatur-Tiefenanalyse",
  "schedule": "0 9 * * *",
  "sources": ["pubmed", "arxiv", "biorxiv"],
  "date_window_days": 1,
  "max_papers_per_source": 50,
  "analysis_tier": "S",
  "dedup": true,
  "resumable": true
}
```

### 8.2 Überwachung ausführen

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json --root ./my-workspace
```

Monitor-Eigenschaften:
- **Fortsetzbar**: Nach Unterbrechung an der letzten Position weitermachen
- **Batch-Verarbeitung**: Kleine Mengen pro Durchlauf, vermeidet Speicherüberlauf
- **Idempotent**: Wiederholte Ausführung erzeugt keine Duplikate (ID-basierte Deduplizierung)
- **Fehlerisoliert**: Fehlschlag einer Einzelanalyse beeinträchtigt nicht den Gesamtprozess

### 8.3 Agent-Cron-Job

Bei Verwendung in einer Hermes-Agent-Umgebung kann ein Cron-Job erstellt werden:

```json
{
  "name": "Tägliche Literatur-Tiefenanalyse",
  "schedule": "0 9 * * *",
  "skills": ["literature-learning-suite"],
  "enabled_toolsets": ["terminal", "file", "web"],
  "workdir": "/pfad/zur/literature-learning-suite"
}
```

**Wichtig**: In unbeaufsichtigten Cron-Läufen keine interaktiven Werkzeuge verwenden. Die
bioRxiv-CDP-Extraktion erfordert manuelle Benutzereingriffe und ist für unbeaufsichtigte
Ausführung nicht geeignet.

### 8.4 Zeitzonenkonfiguration

Bei ungenauen Ausführungszeitpunkten die `timezone`-Einstellung in der Agent-Konfiguration
prüfen. Ohne explizite Zeitzone verwendet der Scheduler möglicherweise UTC.

---

## 9. Fehlerbibliothek & Problembehandlung

Nachfolgend sind in der Produktion tatsächlich aufgetretene Fehler, nach Schweregrad geordnet,
dokumentiert.

### TRAP-001: .pyc-Bytecode-Geist
- **Symptom**: Quellcode-Änderungen zeigen keine Wirkung
- **Ursache**: Python lädt alte .pyc-Dateien aus `__pycache__/`
- **Prävention**: Nach jeder Änderung an `gen_edges.py` ausführen:
  `rm -rf scripts/__pycache__ && python -B`
- **Häufigster Kontext**: Änderungen am Strategiecode in `gen_edges.py`

### TRAP-002: Sandbox-Schreibisolation
- **Symptom**: In `execute_code` geschriebene Dateien erscheinen nicht auf der Festplatte
- **Ursache**: Manche Agent-Umgebungen verwenden temporäre Sandboxen, die Dateien nach
  Abschluss verwerfen
- **Prävention**: `write_file`-Werkzeug oder `terminal` + Python für persistente Schreibvorgänge
- **Muster**: `write_file(pfad, inhalt)` → `terminal("python pfad/zum/skript.py")`

### TRAP-003: .db-Datei wird abgelehnt
- **Symptom**: `read_file` meldet `Cannot read binary file 'papers.db'`
- **Ursache**: Die Erweiterung `.db` wird als Binärdatei erkannt, obwohl es NDJSON-Text ist
- **Prävention**: .db-Dateien immer über `terminal` + `python -c` oder `head` lesen

### TRAP-004: Unicode SyntaxError
- **Symptom**: `SyntaxError: invalid character '→' (U+2192)`
- **Ursache**: Pfeile (→), typografische Anführungszeichen (""), griechische Buchstaben (ΔΨ)
  werden in manchen Shells fehlinterpretiert
- **Prävention**: Nicht-ASCII-Zeichen in Python-Skripten vermeiden. Ersetzungsregeln:
  - `→` → `->`
  - `""` → `'`
  - `ΔΨ` → `Delta`/`Psi`

### TRAP-005: arXiv API stummer Fehlschlag
- **Symptom**: 0 Ergebnisse, keine Fehlermeldung
- **Ursache**: HTTPS-Zugriff auf `export.arxiv.org` liefert 301-Weiterleitung; ohne `-L`
  stiller Fehlschlag
- **Prävention**: `curl -sL "http://export.arxiv.org/api/query?..."` (HTTP + `-L`)

### TRAP-006: PubMed-Proxy-Verlangsamung
- **Symptom**: PubMed-Abfragen laufen in Timeout oder sind extrem langsam
- **Ursache**: PubMed-API ist direkt meist schneller (~780ms); Proxys können verlangsamen
- **Prävention**: Vor PubMed-Abfragen: `unset http_proxy https_proxy`

### TRAP-007: Leerhüllen-S-Publikationen
- **Symptom**: `analysis_tier: "S"` aber T2–T7-Felder sind leer
- **Ursache**: Nur das Label geändert, keinen substanziellen Inhalt geschrieben
- **Erkennung**: `selfcheck_knowledge_graph.py` ausführen, `s_empty`-Zähler prüfen
- **Prävention**: Nach jedem Schreibvorgang sofort das Post-Write-Verifikationsskript
  ausführen (siehe Abschnitt 4.9)

### TRAP-008: Nicht-biologische Kantenverschmutzung
- **Symptom**: `gen_edges.py` produziert `same_journal`-Kanten
- **Ursache**: Vom LLM geschriebene T7-cross_refs können nicht-biologische Beziehungen enthalten
- **Prävention**: Strategie 1 hat eingebauten `NON_BIO_RELS`-Filter
- **Verifikation**: Kanten-Audit-Dimension in `selfcheck_knowledge_graph.py`

### TRAP-009: bioRxiv JATS XML 403
- **Symptom**: `curl https://www.biorxiv.org/content/10.1101/XXXX.source.xml` liefert 403
- **Ursache**: bioRxiv/medRxiv blockiert programmatischen Zugriff auf Quell-XML- und
  PDF-Endpunkte (Stand 2026)
- **Workaround**: Chrome-CDP-Extraktion (`extract_biorxiv_cdp.mjs`) oder API-Abstract
  (300–500 Wörter)

### TRAP-010: Journal-IF-Fuzzy-Match-Fehler
- **Symptom**: IF von „Cell" wird fälschlich als IF von „Cell Reports" zugeordnet
- **Ursache**: Teilstring-Abgleich mit zu grober Granularität
- **Behebung**: Längenverhältnis-Schwellenwert 0,7 in `kg_core.py`
  (`ratio = min(len(norm),len(key)) / max(len(norm),len(key))`)

---

## 10. Plattformspezifische Hinweise

### Linux / macOS

```bash
# Befehlspräfix
python3 scripts/init_workspace.py

# Chrome CDP starten
google-chrome --remote-debugging-port=9223
# oder
chromium --remote-debugging-port=9223

# Gen-Wörterbuch-Regeneration (benötigt R + Bioconductor)
Rscript scripts/export_bioc_genes.R ./my-workspace
```

### Windows

```bash
# Empfohlen: git-bash (MSYS) oder WSL
# Befehlspräfix
python scripts/init_workspace.py   # nicht python3

# Chrome CDP starten
./scripts/biorxiv_chrome_cdp_launcher.bat

# Pfadformate
/c/Users/.../my-workspace   # MSYS-Stil
C:\path\to\my-workspace   # Windows-nativer Stil
# Beide werden akzeptiert
```

**Hinweis**: Bei der Ausführung von `python -c "..."`-Inline-Code in Windows-MSYS-bash
nicht-ASCII-Zeichen in Zeichenketten vermeiden (siehe TRAP-004).

### Plattformübergreifend

- Gebündelte Gen-/Pathway-Wörterbücher und JCR-Journal-Metriken sind sofort einsatzbereit
- `.db`-Dateien sind reiner NDJSON-Text, mit jedem Editor zu öffnen
- Alle Skripte verwenden `pathlib.Path` und behandeln Pfadtrennzeichen automatisch

---

## 11. Bewährte Praktiken

### 11.1 Suchstrategie

- **Erst eng, dann weit**: Mit präzisen Begriffen beginnen, dann schrittweise lockern
- **Erst PubMed, dann arXiv**: PubMed antwortet schneller (~780ms direkt) und hat breitere
  Abdeckung
- **Jede Suche protokollieren**: Vollständige Abfragezeichenkette, Ausführungszeit,
  Trefferanzahl — essenziell für Reproduzierbarkeit und Audit
- **Null Ergebnisse ≠ Nichtexistenz**: Zuerst Syntax, Konnektivität, Datumsfilter und
  Ratenbegrenzungen prüfen

### 11.2 Analysestrategie

- **Volltext vor Abstract**: Volltext → vollständige 7-Schichten-Analyse. Nur Abstract →
  Aussagekraft einschränken, `abstract_only` kennzeichnen
- **Kein Batch**: Jede Publikation einzeln analysieren. Vorlagen-Kopieren ist die
  Hauptursache für Leerhüllen-S
- **Evidenz muss konkret sein**: „2847 NSCLC, SHAP, EFS HR=0,52" nicht „Autoren bewiesen"
- **Mechanismus präzise bis zur Modifikationsstelle**: „NF-κB p65 Ser536-Phosphorylierung"
  nicht „aktiviert Signalweg"

### 11.3 Wartungsstrategie

- **Täglich `selfcheck_knowledge_graph.py`**: Als letzten Schritt des Cron-Jobs, vor der
  Berichterstellung auditieren
- **Regelmäßig .pyc bereinigen**: Besonders nach Änderungen an `gen_edges.py`
- **Wöchentliches Backup**: papers.db, edges.db, concepts.db sind die drei kritischsten
  Dateien
- **Protokolle überwachen**: `logs/`-Verzeichnis auf Anomalien prüfen

### 11.4 Teamarbeit

- **Ein Arbeitsbereich pro Person**: Jeder Forschende pflegt seinen eigenen Arbeitsbereich,
  um Nebenläufigkeitskonflikte zu vermeiden
- **Gemeinsame Konzeptbibliothek**: concepts.db kann arbeitsbereichübergreifend
  zusammengeführt werden
- **Kanten sind reproduzierbar**: edges.db wird vollständig von gen_edges.py generiert
  und kann jederzeit neu erstellt werden
- **Git verwenden**: .db-Dateien (NDJSON-Text) können versioniert werden;
  `git diff` zeigt Änderungen

---

## 12. FAQ

**F: Warum NDJSON statt SQLite?**

A: NDJSON ist menschenlesbar, versionskontrollfreundlich, CLI-abfragbar und hat keine externen
Abhängigkeiten. Für „Einmal schreiben, oft lesen, gelegentlich anhängen"-Workloads ist NDJSON
die passendere Wahl.

**F: Wie handhabe ich inkrementelle Aktualisierungen?**

A: `add_paper()` und `kg.py add` haben eingebaute Deduplizierung anhand stabiler IDs. Einfach
neue Datensätze anhängen; bestehende werden nicht überschrieben.

**F: Wie oft soll gen_edges.py ausgeführt werden?**

A: Nach jeder Charge neuer Publikationen. Im täglichen Cron-Job als letzten Schritt ausführen
(.pyc bereinigen → gen_edges → build_network → selfcheck).

**F: Kann ich eigene Journal-IF-Daten verwenden?**

A: Ja. Ihre IF-Daten als NDJSON in `my-workspace/journal_metrics.db` ablegen (gleiche Felder
wie `assets/data/journal_metrics_2024.json`). Das System bevorzugt Workspace-Daten und fällt
auf die gebündelten JCR-2024-Daten zurück.

**F: Wie beschaffe ich bioRxiv-Volltexte?**

A: API-Abstracts (300–500 Wörter) sind immer verfügbar. Für JavaScript-gerenderte Volltexte
wird die Chrome-CDP-Extraktion benötigt (`extract_biorxiv_cdp.mjs`). Dieser Pfad erfordert
manuelle Browser-Bedienung zur Cloudflare-Verifikation.

**F: Werden deutschsprachige Publikationen unterstützt?**

A: Die aktuellen Suchquellen (PubMed/arXiv/bioRxiv/Crossref) sind überwiegend
englischsprachig. Deutschsprachige Publikationen können manuell erstellt und via
`normalize_records.py` + `kg.py add` eingepflegt werden. Das Analyseprotokoll selbst ist
sprachunabhängig.

**F: Wie aktualisiere ich ältere Analyse-Datensätze?**

A: Siehe `references/s-tier-upgrade-workflow.md`. Ablauf: Auditieren → Batch-Upgrade
(10 pro Satz) → Kanten und Netzwerk neu aufbauen → Erneut auditieren.

---

## 13. Glossar

| Deutsch | Englisch | Erläuterung |
|---------|----------|-------------|
| Wissensgraph | Knowledge Graph | NDJSON-basiertes Netzwerk aus Publikationen, Konzepten und semantischen Kanten |
| S-Tier-Analyse | S-tier Analysis | Vollständiges 7-Schichten-Dissektionsprotokoll |
| Leerhülle-S | Empty-shell S | Als S markierter Datensatz ohne substanziellen T2–T7-Inhalt |
| Behauptung-Evidenz-Synthese-Kette (BES) | Claim-Evidence-Synthesis Chain (CES) | Kerneinheit der T3-Analyse |
| Verborgene Organisationsachse | Hidden Organizing Axis | T5-Schicht — Muster, die die Publikation nicht explizit benennt |
| Semantische Kante | Semantic Edge | Biologisch bedeutsame Verbindung zwischen Publikationen |
| NDJSON | Newline-Delimited JSON | Speicherformat mit einem vollständigen JSON-Objekt pro Zeile |
| Invertierter Index | Inverted Index | Kerndatenstruktur von gen_edges.py (O(M)-Komplexität) |
| Arbeitsbereich | Workspace | Laufzeit-Datenverzeichnis (my-workspace/) |
| Stabiler Identifikator | Stable Identifier | Persistent zitierbare ID (PMID/DOI/arXiv-ID) |
| Evidenzbewertungsmatrix | Evidence Rubric | Standardisiertes Bewertungssystem basierend auf Verzerrungsrisiko/Stichprobengröße/Reproduzierbarkeit |
| JCR | Journal Citation Reports | Clarivate-Journalzitationsberichte |
| CDP | Chrome DevTools Protocol | Browser-Automatisierungsprotokoll |
| Zweistufiger Abgleich | Two-stage Matching | Regex-Kandidatenextraktion → Set-Lookup-Bestätigung |

---

> Verwaltet vom Literature Learning Suite Projekt.
> Version 1.3.0, letzte Aktualisierung 2026-06-09.
> Lizenz: CC BY-NC-SA 4.0 (Namensnennung — Nicht-kommerziell — Weitergabe unter gleichen Bedingungen)
