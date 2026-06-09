# Literature Learning Suite — Guía Completa del Usuario

> **Versión**: 1.3.0 | **Licencia**: CC BY-NC-SA 4.0

---

## Índice

1. [Descripción General y Filosofía de Diseño](#1-descripción-general-y-filosofía-de-diseño)
2. [Arquitectura del Sistema](#2-arquitectura-del-sistema)
3. [Inicio Rápido](#3-inicio-rápido)
4. [Metodología Central: Protocolo de 7 Capas de Nivel S](#4-metodología-central-protocolo-de-7-capas-de-nivel-s)
5. [Modelo de Datos del Grafo de Conocimiento](#5-modelo-de-datos-del-grafo-de-conocimiento)
6. [Motor de Generación de Enlaces v3.1](#6-motor-de-generación-de-enlaces-v31)
7. [Referencia de Herramientas](#7-referencia-de-herramientas)
8. [Monitoreo Automatizado](#8-monitoreo-automatizado)
9. [Biblioteca de Trampas y Solución de Problemas](#9-biblioteca-de-trampas-y-solución-de-problemas)
10. [Notas de Plataforma](#10-notas-de-plataforma)
11. [Mejores Prácticas](#11-mejores-prácticas)
12. [Preguntas Frecuentes](#12-preguntas-frecuentes)
13. [Glosario](#13-glosario)

---

## 1. Descripción General y Filosofía de Diseño

### 1.1 Qué es Esto

El Literature Learning Suite (LLS) es un sistema completo para el descubrimiento de literatura académica, el análisis profundo, la construcción de grafos de conocimiento y el monitoreo automatizado. No es "otro gestor bibliográfico más" — es un **sistema operativo de cognición literaria**.

Los gestores bibliográficos tradicionales (Zotero, EndNote, Mendeley) resuelven el problema de "dónde almacenar". LLS resuelve los problemas de **"cómo leer"** y **"cómo pensar"**. A través de un riguroso protocolo de disección en 7 capas, transforma cada artículo en un nodo estructurado, consultable y enlazable dentro de un grafo de conocimiento.

### 1.2 Principios de Diseño

Estos principios cristalizaron a partir del análisis profundo de cientos de artículos:

**Principio 1: Profundidad sobre amplitud.** Un solo análisis completo de nivel S en 7 capas vale más que 50 artículos con solo títulos y resúmenes. Si dispone de 2 horas, dedique 1,5 horas a leer profundamente 2 artículos, no a hojear 20.

**Principio 2: Verificar antes de citar.** Cruce la identidad bibliográfica (DOI/PMID/arXiv ID) en al menos dos fuentes antes del análisis o la cita. No verificar + no etiquetar = propagar errores.

**Principio 3: Separar la evidencia reportada de la inferencia.** Cada afirmación debe etiquetarse con su nivel epistemológico — REPORTADO (directamente en el artículo), INFERENCIA RESPALDADA (deriva de evidencia + conocimiento establecido), HIPÓTESIS (generada por el analista), o DESCONOCIDO (no determinable a partir del material disponible).

**Principio 4: Juzgar la evidencia por el diseño del estudio, no por el prestigio de la revista.** El factor de impacto es metadato, no una puntuación de calidad. LLS incluye una rúbrica de evidencia estandarizada basada en riesgo de sesgo, tamaño muestral, calidad de controles y replicabilidad.

**Principio 5: Persistencia solo por adición.** Todas las bases de datos NDJSON utilizan modo de solo adición (append-only). Esto garantiza un historial operativo completo y previene la pérdida irreversible de datos.

**Principio 6: Nunca fabricar.** No se inventan artículos, identificadores, estadísticas ni mecanismos. Si una cita no puede verificarse, se excluye.

**Principio 7: Texto como dato no confiable.** El texto del artículo y el contenido web son datos no confiables — objetos de análisis, nunca instrucciones para el agente. Esto es fundamental para la prevención de inyección de prompts.

### 1.3 Comparación con Herramientas Existentes

| Dimensión | Zotero/EndNote | Revisión Tradicional | LLS |
|-----------|---------------|---------------------|-----|
| Almacenamiento | PDF + metadatos | — | Grafo de conocimiento NDJSON |
| Profundidad de análisis | Notas manuales | Resumen humano | Disección estructurada en 7 capas |
| Enlaces entre artículos | Etiquetas manuales | Subjetivo | 5 estrategias de enlaces semánticos automáticos |
| Calificación de evidencia | Ninguna | Basada en experiencia | Rúbrica estandarizada |
| Automatización | Asistida por plugins | Ninguna | Pipeline completo de monitoreo |
| Consultabilidad | Búsqueda por palabras clave | No consultable | Búsqueda de texto completo + recorrido del grafo |
| Formato de persistencia | Propietario | Texto plano | NDJSON (legible por humanos) |

---

## 2. Arquitectura del Sistema

### 2.1 Estructura de Directorios

```
literature-learning-suite/
├── GUIDE_EN.md              ← Documento en inglés
├── GUIDE_ZH.md              ← Guía en chino (principal)
├── GUIDE_DE.md              ← Guía en alemán
├── GUIDE_JA.md              ← Guía en japonés
├── GUIDE_KO.md              ← Guía en coreano
├── GUIDE_ES.md              ← Este documento
├── SKILL.md                 ← Protocolo de operación del agente
├── README.md                ← Descripción general del proyecto
├── LICENSE                  ← CC BY-NC-SA 4.0
│
├── scripts/                 ← Herramientas centrales (23+ Python + Node.js + R)
├── assets/data/             ← Datos de referencia validados (incluidos)
│   ├── bioc_genes.json      ← 90.125 símbolos génicos humanos y de ratón
│   ├── kegg_pathways.json   ← 25.939 términos de vías KEGG + GO BP
│   └── journal_metrics_2024.json ← 21.800 revistas con FI/cuartiles
├── references/              ← Documentación de metodología y protocolos (25+ docs)
└── tests/                   ← Pruebas unitarias + datos sintéticos
```

### 2.2 Pipeline de Flujo de Datos

```
Pregunta de Investigación
       │
       ▼
  [1. Formular estrategia de búsqueda]  ← Marco PICO/PECO
       │
       ▼
  [2. Búsqueda multi-fuente]    ← PubMed/arXiv/bioRxiv/Crossref
       │
       ▼
  [3. Verificar identidad]      ← Verificación cruzada DOI/PMID
       │
       ▼
  [4. Obtener texto completo]   ← PMC XML / arXiv HTML / extracción CDP
       │
       ▼
  [5. Análisis de 7 capas Nivel S] ← Razonamiento profundo LLM (T1-T7)
       │
       ▼
  [6. Calificar evidencia]      ← Riesgo de sesgo / tamaño muestral / replicabilidad
       │
       ▼
  [7. Persistir grafo de conocimiento] ← NDJSON append
       │
  ┌────┼────┐
  ▼    ▼    ▼
 [Enlaces] [Resumen] [Red]
       │
       ▼
  [9. Autoverificación de calidad + Monitoreo]
```

### 2.3 Espacio de Trabajo en Tiempo de Ejecución

Creado por `init_workspace.py`, un directorio autónomo:

```
my-workspace/
├── workspace.json           ← Metadatos del espacio de trabajo
├── papers.db                ← NDJSON, un artículo por línea
├── concepts.db              ← NDJSON, un concepto por línea
├── edges.db                 ← NDJSON, un enlace semántico por línea
├── queries.db               ← Registro de búsquedas
├── data/                    ← Diccionarios de genes/vías sembrados
├── fulltext/                ← Documentos descargados
├── daily_digest/            ← Resúmenes diarios generados
├── config/                  ← Configuración del monitor
├── exports/                 ← Exportaciones BibTeX
└── logs/                    ← Registros de ejecución
```

---

## 3. Inicio Rápido

### 3.1 Requisitos

- Python 3.10+
- pip
- Node.js v24+ (solo para extracción CDP de bioRxiv)
- R + Bioconductor (solo para regeneración del diccionario génico)

### 3.2 Instalación

```bash
git clone <repo-url>
cd literature-learning-suite
pip install -r scripts/requirements.txt
```

### 3.3 Inicializar el Espacio de Trabajo

```bash
python scripts/init_workspace.py --root ./my-workspace
```

Esto crea la estructura de directorios del espacio de trabajo y siembra el diccionario génico (90.125 símbolos) y las etiquetas de vías (25.939 términos).

### 3.4 Establecer Variable de Entorno

```bash
# Linux/macOS
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"

# Windows PowerShell
$env:LITERATURE_KG_ROOT = (Resolve-Path ./my-workspace)
```

Si no se define, los scripts usan por defecto `./literature-workspace`.

### 3.5 Primera Búsqueda

```bash
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10
python scripts/literature_search.py arxiv "single cell foundation model" -n 10
python scripts/literature_search.py biorxiv "spatial transcriptomics" -n 10
```

### 3.6 Ingesta de Artículos

```bash
python scripts/normalize_records.py search-results.jsonl -o normalized.jsonl
python scripts/validate_records.py normalized.jsonl
python scripts/kg.py --root ./my-workspace add normalized.jsonl
python scripts/kg.py --root ./my-workspace stats
```

---

## 4. Metodología Central: Protocolo de 7 Capas de Nivel S

Este es el corazón de LLS. Cada artículo de investigación debe someterse al análisis completo de 7 capas.

### 4.1 Por Qué 7 Capas

La lectura tradicional se detiene en "leer resumen → etiquetar → escribir dos frases". Este procesamiento superficial impide capturar la estructura lógica profunda, descubrir conexiones implícitas entre artículos, evaluar la verdadera solidez de la evidencia o formar conocimiento estructurado consultable.

El protocolo de 7 capas descompone la lectura en siete dimensiones cognitivas independientes pero complementarias, obligando al analista (humano o LLM) a producir resultados sustantivos en cada una.

### 4.2 Capa 1 (T1): Perfil Bibliográfico y del Estudio

Establecer identidad verificable y contexto: título completo, primeros autores y autor de correspondencia con afiliaciones, revista (anotada automáticamente con el FI de JCR 2024 y cuartil mediante `enrich_paper_if()`), identificadores estables (DOI/PMID/arXiv ID), tipo de documento, diseño del estudio, tamaño muestral, estado del texto completo.

Esta capa puede automatizarse.

### 4.3 Capa 2 (T2): Pregunta Científica Central

Destilar el artículo en **una pregunta central falsable** y descomponerla en **≥5 subpreguntas comprobables**. La pregunta central debe ser mecanicista, no descriptiva.

**Ejemplo (bueno):**
> Pregunta central: ¿Cómo induce GDF15 la disfunción de células NK a través de la señalización de receptores xenobióticos?
> Subpreguntas: ¿Qué tumores expresan GDF15? ¿Qué receptores median la señal? ¿Qué cascada de señalización deteriora la citotoxicidad NK? ¿Es reversible la disfunción? ¿El bloqueo de GDF15 restaura la inmunidad antitumoral in vivo?

### 4.4 Capa 3 (T3): Cadenas Afirmación-Evidencia-Síntesis (≥5)

Transformar los hallazgos centrales en ≥5 cadenas C-E-S independientes con evidencia concreta.

| Campo | Requisito | Umbral |
|-------|-----------|--------|
| `claim` (afirmación) | Conclusión falsable en sus propias palabras | — |
| `evidence` (evidencia) | Datos concretos: N, tamaño del efecto, valor p, ensayo | **>20 caracteres** (disparador de caparazón vacío) |
| `synthesis` (síntesis) | Cómo la evidencia respalda/debilita la afirmación; vínculos interdisciplinarios | — |
| `strength` (solidez) | 1-5 estrellas, justificado por calidad de evidencia | — |
| `uncertain` (incertidumbre) | Explicaciones alternativas, validaciones faltantes | — |

**Bueno:** "2847 pacientes con NSCLC, análisis de interacción SHAP, SLE HR=0,52 (IC 95% 0,34-0,78), p=0,002, revisión radiológica independiente."

**Malo (caparazón vacío):** "Los autores demostraron su hipótesis."

### 4.5 Capa 4 (T4): Cascada de Mecanismos

Mapear la cadena causal completa desde el desencadenante hasta el fenotipo, con precisión hasta los sitios de modificación:

```
Desencadenante → Receptor → Segundo mensajero/quinasa → Factor de transcripción → Gen diana → Fenotipo
```

Debe incluir: ≥3 pasos en la cascada, sitios clave de modificación (ej. "fosforilación de NF-κB p65 Ser536", no "activación de NF-κB"), efectos posteriores, ≥1 bucle de retroalimentación. Etiquetar cada paso: `demostrado`, `respaldado`, `conocimiento previo`, o `hipótesis`.

### 4.6 Capa 5 (T5): Ejes Organizativos Ocultos (≥3)

Descubrir patrones que el artículo **no** declara explícitamente — supuestos implícitos, organización espaciotemporal, efectos de selección o lógica experimental más profunda.

Cada eje: una `observación` (hecho verificable del artículo) + una `interpretación` (el patrón más profundo descubierto por el analista).

T5 evalúa su capacidad de síntesis, no su habilidad para copiar la sección de Discusión.

**Ejemplo:**
> Observación 1: En todos los análisis de reducción dimensional, las muestras del borde tumoral se separan consistentemente de las del núcleo como un clúster independiente.
> Interpretación 1: El estudio define implícitamente el "borde" y no el "núcleo" como el compartimento determinante de la enfermedad, lo que explica por qué las firmas génicas del núcleo tienen paradójicamente menor valor pronóstico — la heterogeneidad molecular del núcleo queda enmascarada por las señales microambientales del borde.
>
> Observación 2: En el análisis unicelular, los cambios en las subpoblaciones inmunitarias preceden a los cambios en las subpoblaciones tumorales.
> Interpretación 2: Los datos de series temporales del estudio sugieren un modelo de progresión de la enfermedad "primero inmunitario" — la remodelación del microambiente es el motor de la evolución tumoral, no su consecuencia. Si esta interpretación es correcta, los puntos de intervención temprana deberían dirigirse a las células inmunitarias, no a las tumorales.

### 4.7 Capa 6 (T6): Contribución Conceptual

Identificar: nuevos conceptos con definiciones operativas (≥1), visiones previas desafiadas o acotadas (≥1), avances metodológicos (≥1), y condiciones de contorno (¿cuándo NO se generaliza esto?).

**Ejemplo:**
> Nuevo concepto:
> - "Punto de control inmunometabólico" (Immunometabolic Checkpoint): definido como un eje de inmunosupresión mediado por metabolitos a través de receptores inmunitarios, distinto de los puntos de control clásicos basados en ligandos proteicos.
>   - Definición operativa: debe cumplir simultáneamente (a) metabolito pequeño como ligando, (b) señalización a través de receptor inmunitario, (c) supresión funcional reversible de células inmunitarias.
>
> Visión previa desafiada:
> - Corrige la hipótesis de función única de "GDF15 como factor anorexigénico", demostrando que posee un rol independiente en la inmunidad tumoral.
>
> Avance metodológico:
> - Establece un método de detección de unión metabolito-receptor independiente de anticuerpos (marcaje SILAC + entrecruzamiento químico + espectrometría de masas).
>
> Condiciones de contorno:
> - El mecanismo solo opera en tumores con alta expresión de GDF15 (≥ 2 veces la mediana).

### 4.8 Capa 7 (T7): Relaciones entre Artículos (≥5)

Conectar el artículo con otros registros verificados usando relaciones biológicas sustantivas (`supports`, `contradicts`, `extends`, `replicates`, `shared_mechanism`, etc.). Cada relación requiere una descripción de 60-150 palabras explicando POR QUÉ.

**Explícitamente prohibidos:** `same_journal`, `same_issue`, `same_author` (ruido no biológico).

### 4.9 Detección de Caparazón Vacío (Empty-Shell S)

Un registro etiquetado `S` es un caparazón vacío si se cumple CUALQUIERA de:
- `tier2_subquestions` vacío
- `tier3_ces_chains` < 5
- Algún campo `evidence` ≤ 20 caracteres
- `tier4` cascada < 3 pasos
- `tier5_hidden_axis` < 3
- `tier7_cross_refs` < 5
- Alguna relación T7 está en la lista de prohibidas

**Verificación post-escritura:**
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

### 4.10 Artículos No Científicos

Artículos clasificados como noticias/editoriales/erratas/retractaciones: etiquetar `analysis_tier: "NR"`, registrar el tipo real en `core_findings`, y NO fabricar contenido de análisis.

---

## 5. Modelo de Datos del Grafo de Conocimiento

### 5.1 Justificación del Diseño

LLS utiliza NDJSON (JSON delimitado por saltos de línea) en lugar de SQLite porque:
1. Legible y editable por humanos
2. Compatible con control de versiones (git diff funciona línea por línea)
3. Solo adición para integridad de datos
4. Consultable con herramientas CLI estándar (head, grep, jq)
5. Cero dependencias externas de bases de datos

### 5.2 papers.db

Un registro completo de artículo por línea. Campos obligatorios: `id`, `title`, `source`, `retrieved_at`. El esquema completo S-tier v4.0 incluye campos de análisis T2-T7 más `entities` para anotaciones de genes/vías/tipos celulares/enfermedades.

### 5.3 concepts.db

```json
{"id": "CONCEPT:punto_control_inmunometabolico", "name": "Punto de Control Inmunometabólico",
 "type": "mechanism", "definition": "...", "source_papers": ["PMID:42251595"]}
```

### 5.4 edges.db

```json
{"source": "PMID:42251595", "target": "PMID:39988000", "relation": "extends",
 "description": "Justificación biológica de 60-150 palabras", "provenance": "analyst"}
```

### 5.5 queries.db

Registro de auditoría de búsquedas: `source`, `query`, `executed_at`, `result_count`, `parameters`.

---

## 6. Motor de Generación de Enlaces v3.1

### 6.1 Descripción General

`gen_edges.py` construye conexiones semánticamente significativas entre artículos utilizando 5 estrategias complementarias, todas implementadas con índices invertidos (complejidad O(M), no O(N²)).

### 6.2 Cinco Estrategias

| # | Estrategia | Relación | Lógica |
|---|-----------|----------|--------|
| 1 | Refs explícitas | extends, accompanied_by | tier7_cross_refs del análisis profundo (filtrado no-bio) |
| 2 | Moléculas compartidas | shares_molecules (≥2) | 90.125 genes, emparejamiento en dos etapas |
| 2.5 | Solapamiento textual | shares_topic (≥4) | Bolsa de palabras de hallazgos centrales → índice invertido |
| 3 | Enfermedad × método | shares_disease_method | Producto cruzado: etiqueta enfermedad × etiqueta método |
| 4 | Eje oculto | shares_paradigm | Resonancia de patrones profundos de T5 |
| 5 | Nodos conceptuales | defines_concept | Artículo → enlaces a conceptos desde concepts.db |

### 6.3 Uso

```bash
# CRÍTICO: limpiar caché de bytecode antes de cada ejecución
rm -rf scripts/__pycache__
python -B scripts/gen_edges.py
python scripts/build_network.py  # actualizar visualización
```

---

## 7. Referencia de Herramientas

### Búsqueda y Descubrimiento

- `literature_search.py` — Búsqueda multi-fuente (pubmed/arxiv/crossref/biorxiv)
- `search_arxiv.py` — Búsqueda específica de arXiv con filtros por autor/categoría
- `download_biorxiv_api.py` — Descarga por lotes de bioRxiv/medRxiv

### Verificación y Normalización

- `verify_citation.py` — Verificación cruzada de identidad DOI/PMID/arXiv
- `normalize_records.py` — Estandarización + desduplicación
- `validate_records.py` — Validación de esquema JSON

### Adquisición de Texto Completo

- `fulltext_fetch.py` — Descargador unificado (PMC XML + arXiv HTML)
- `extract_pymupdf.py` — Extracción de texto y tablas de PDF
- `extract_marker.py` — Extracción OCR de PDF
- `extract_biorxiv_cdp.mjs` — Extracción CDP de Chrome para bioRxiv/medRxiv

### Grafo de Conocimiento

- `kg.py` — CLI: add/stats/search/audit
- `kg_core.py` — Biblioteca: CRUD + anotación automática de FI mediante `enrich_paper_if()`

### Análisis y Salida

- `gen_edges.py` — Generador de enlaces semánticos v3.1
- `gen_digest.py` — Resumen diario en Markdown
- `build_network.py` — Grafo HTML interactivo de fuerza dirigida
- `selfcheck_knowledge_graph.py` — Auditoría de calidad en 10 dimensiones
- `export_citations.py` — Exportación BibTeX/CSL

### Automatización

- `monitor.py` — Monitor por lotes reanudable
- `init_workspace.py` — Inicializar un nuevo espacio de trabajo

---

## 8. Monitoreo Automatizado

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json --root ./my-workspace
```

Características: reanudable, procesamiento por lotes, idempotente (desduplicación por ID estable), tolerante a fallos aislados.

Para cron en agente Hermes:
```json
{"name": "Análisis Profundo Diario de Literatura", "schedule": "0 9 * * *",
 "skills": ["literature-learning-suite"], "enabled_toolsets": ["terminal", "file", "web"]}
```

⚠️ Evitar herramientas interactivas en ejecuciones cron desatendidas. La extracción CDP de bioRxiv requiere operación manual.

---

## 9. Biblioteca de Trampas y Solución de Problemas

### TRAP-001: Fantasma de Bytecode .pyc
**Síntoma:** Los cambios en el código fuente no surten efecto.
**Prevención:** `rm -rf scripts/__pycache__ && python -B`

### TRAP-002: Aislamiento de Escritura en Sandbox
**Síntoma:** Archivos escritos en `execute_code` no aparecen en disco.
**Prevención:** Usar la herramienta `write_file` o `terminal` + Python.

### TRAP-003: Archivo .db Rechazado por read_file
**Síntoma:** `Cannot read binary file 'papers.db'`
**Prevención:** Usar `terminal` + `python -c` o `head` para leer archivos .db.

### TRAP-004: SyntaxError por Caracteres Unicode
**Síntoma:** `SyntaxError: invalid character '→'`
**Prevención:** Reemplazar `→` por `->`, `""` por `'`, `ΔΨ` por `Delta`/`Psi`.

### TRAP-005: Falla Silenciosa de arXiv
**Síntoma:** Cero resultados, sin error.
**Prevención:** Usar `curl -sL "http://export.arxiv.org/api/query?..."` (HTTP + -L).

### TRAP-006: Ralentización por Proxy en PubMed
**Síntoma:** Las consultas a PubMed expiran por tiempo.
**Prevención:** `unset http_proxy https_proxy` antes de llamadas a PubMed.

### TRAP-007: Registros S de Caparazón Vacío
**Síntoma:** `analysis_tier: "S"` pero T2-T7 vacíos.
**Detección:** Ejecutar `selfcheck_knowledge_graph.py`, verificar el conteo `s_empty`.

### TRAP-008: Contaminación por Enlaces No Biológicos
**Síntoma:** `gen_edges.py` produce enlaces `same_journal`.
**Prevención:** La estrategia 1 tiene un filtro integrado `NON_BIO_RELS`.

### TRAP-009: Error 403 en JATS XML de bioRxiv
**Síntoma:** El endpoint `.source.xml` devuelve 403.
**Solución alternativa:** Extracción CDP de Chrome o resúmenes API (300-500 palabras).

### TRAP-010: Error de Coincidencia Difusa del FI de Revista
**Síntoma:** El FI de "Cell" se asigna al de "Cell Reports".
**Solución:** Umbral de relación de longitud > 0,7 en `kg_core.py`.

---

## 10. Notas de Plataforma

### Linux / macOS
- `python3` es el comando habitual
- Chrome CDP: `google-chrome --remote-debugging-port=9223`
- Regeneración de genes: `Rscript scripts/export_bioc_genes.R`

### Windows
- Usar git-bash (MSYS) o WSL
- `python` (no `python3`)
- Chrome CDP: `scripts/biorxiv_chrome_cdp_launcher.bat`
- Evitar caracteres no ASCII en instrucciones Python de una línea ejecutadas en MSYS (ver TRAP-004)

### Todas las Plataformas
- Los diccionarios de genes/vías incluidos y las métricas JCR 2024 funcionan sin configuración adicional
- `enrich_paper_if()` anota automáticamente el FI desde la base de datos de 21.800 revistas
- Los archivos `.db` son texto NDJSON plano

---

## 11. Mejores Prácticas

### Estrategia de Búsqueda
- Comenzar con búsquedas precisas y luego ampliar
- PubMed primero (más rápido, ~780ms directo), luego arXiv
- Registrar cada búsqueda para reproducibilidad
- Cero resultados ≠ evidencia de ausencia — verificar sintaxis, conectividad, filtros

### Estrategia de Análisis
- Texto completo > resumen. Solo resumen → limitar afirmaciones, etiquetar `abstract_only`
- No procesar por lotes. Un artículo a la vez. Copiar plantillas es la causa principal de caparazones vacíos
- La evidencia debe ser concreta: "2847 NSCLC, SHAP, SLE HR=0,52" no "los autores demostraron"
- Mecanismo preciso hasta el sitio: "fosforilación de NF-κB p65 Ser536" no "activa la vía"

### Mantenimiento
- Ejecutar `selfcheck_knowledge_graph.py` diariamente como último paso del cron
- Limpiar .pyc después de cada cambio en gen_edges.py
- Respaldo semanal de papers.db, edges.db, concepts.db
- Monitorear logs/ en busca de anomalías

### Colaboración en Equipo
- Un espacio de trabajo por investigador para evitar conflictos de escritura concurrente
- concepts.db puede fusionarse entre espacios de trabajo
- edges.db se regenera completamente con gen_edges.py, puede reconstruirse en cualquier momento
- Usar Git: los archivos .db (texto NDJSON) pueden versionarse; git diff muestra los cambios

---

## 12. Preguntas Frecuentes

**P: ¿Por qué NDJSON en lugar de SQLite?**
R: Legible por humanos, compatible con Git, consultable con CLI, cero dependencias. Para cargas de trabajo de "escribir una vez, leer muchas, añadir ocasionalmente", NDJSON es más adecuado.

**P: ¿Cómo gestionar actualizaciones incrementales?**
R: `add_paper()` y `kg.py add` incluyen desduplicación integrada por ID estable. Simplemente añada nuevos registros.

**P: ¿Con qué frecuencia ejecutar gen_edges.py?**
R: Después de cada lote de nuevos artículos. En el cron diario, ejecutar como último paso.

**P: ¿Puedo usar mis propios datos de FI de revistas?**
R: Sí. Coloque sus datos NDJSON en `workspace/journal_metrics.db`. El sistema prefiere los datos del espacio de trabajo y usa los datos JCR 2024 incluidos como respaldo.

**P: ¿Cómo obtener el texto completo de bioRxiv?**
R: Los resúmenes API (300-500 palabras) están siempre disponibles. El texto completo renderizado con JavaScript requiere extracción CDP de Chrome con verificación manual de Cloudflare.

**P: ¿Artículos en español u otros idiomas?**
R: Las fuentes de búsqueda actuales se centran en inglés. Los artículos en otros idiomas pueden crearse manualmente e ingerirse mediante `normalize_records.py` + `kg.py add`. El protocolo de análisis es independiente del idioma.

**P: ¿Cómo actualizar registros de análisis de versiones anteriores?**
R: Consulte `references/s-tier-upgrade-workflow.md`. Flujo: auditar → actualizar por lotes (10/lote) → reconstruir enlaces y red → volver a auditar.

---

## 13. Glosario

| Término | Definición |
|---------|------------|
| Grafo de Conocimiento | Red basada en NDJSON de artículos, conceptos y enlaces semánticos |
| Análisis de Nivel S | Protocolo completo de disección en 7 capas |
| Caparazón Vacío (Empty-shell S) | Registro etiquetado S sin contenido sustantivo en T2-T7 |
| Cadena Afirmación-Evidencia-Síntesis (CES) | Unidad central del análisis T3 |
| Eje Organizativo Oculto | Capa T5 — patrones que el artículo no declara explícitamente |
| Enlace Semántico | Conexión biológicamente significativa entre artículos |
| NDJSON | Formato de almacenamiento JSON delimitado por saltos de línea |
| Índice Invertido | Estructura de datos central de gen_edges.py (complejidad O(M)) |
| Espacio de Trabajo | Directorio de datos en tiempo de ejecución (my-workspace/) |
| Identificador Estable | ID persistentemente citable (PMID/DOI/arXiv ID) |
| Rúbrica de Evidencia | Sistema estandarizado de calificación de evidencia |
| CDP | Protocolo Chrome DevTools |
| Emparejamiento en Dos Etapas | Extracción de candidatos por regex → confirmación por búsqueda en conjunto |
| FI (Factor de Impacto) | Métrica de revista según Journal Citation Reports (JCR) |
| PICO/PECO | Marco para formular preguntas de investigación: Población, Intervención/Exposición, Comparación, Resultado (Outcome) |

---

> Mantenido por el proyecto Literature Learning Suite.
> Versión 1.3.0, última actualización 2026-06-09.
> Licencia: CC BY-NC-SA 4.0 (Atribución-NoComercial-CompartirIgual)
