# Literature Learning Suite

> 从研究问题到知识图谱的完整文献认知操作系统——不是文献管理工具，而是文献**思考**工具。

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.3.0-green.svg)](VERSION)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](scripts/requirements.txt)

---

## 目录

- [1. 项目定位](#1-项目定位)
- [2. 解决什么问题](#2-解决什么问题)
- [3. 核心能力](#3-核心能力)
- [4. 系统架构](#4-系统架构)
- [5. 完整工作流](#5-完整工作流)
- [6. 安装与快速开始](#6-安装与快速开始)
- [7. 核心方法论预览：S 级 7 层解剖](#7-核心方法论预览s-级-7-层解剖)
- [8. 工具链概览](#8-工具链概览)
- [9. 文档导航](#9-文档导航)
- [10. 打包数据](#10-打包数据)
- [11. 许可证](#11-许可证)
- [12. 目录结构](#12-目录结构)

---

## 1. 项目定位

**Literature Learning Suite**（以下简称 LLS）是一套自包含的学术文献深度分析系统。
它的核心任务不是"文献存在哪里"（Zotero、EndNote 已经解决了这个问题），而是：

> **怎么读**——如何从一篇论文中系统性地提取出结构化的、可验证的、可关联的知识
>
> **怎么想**——如何评估证据强度、发现隐藏假设、建立跨论文的语义关联

LLS 将文献阅读从"手动标记 + 写几句总结"升级为**7 层结构化解剖**，
将分散的论文转化为一个**可查询、可遍历、可自检的知识图谱**。

### 与传统工具的差异

| | Zotero / EndNote | 传统文献综述 | LLS |
|---|---|---|---|
| 存储格式 | PDF + 元数据 | 纯文本 | **NDJSON 知识图谱** |
| 分析深度 | 手动笔记 | 人工总结 | **7 层结构化解剖** |
| 跨文献关联 | 手动标签 | 主观判断 | **5 策略自动语义边** |
| 证据评分 | 无 | 依赖经验 | **标准化评分标准** |
| 自动化 | 插件辅助 | 无 | **完整监控流水线** |
| 可查询性 | 关键词搜索 | 不可查询 | **全文检索 + 图遍历** |
| 持久化 | 专有格式 | 脆弱 | **NDJSON（人可读写、git 友好）** |
| 影响因子 | 无 | 手动查 | **21,800 期刊自动标注** |

---

## 2. 解决什么问题

### 问题 1：文献阅读停留在表面

大多数研究者读论文的方式：看标题 → 读摘要 → 标记"有用" → 写两句话。这种浅层处理无法捕捉：

- 论文的**深层逻辑结构**（他们为什么这样设计实验？）
- 论文**没有明说**的假设和限制
- 多条证据链之间的**独立性和互补性**
- 证据的**真实强度**（p 值不等于效果量，统计显著不等于临床显著）

**LLS 的解决方案**：强制性的 7 层解剖协议，每层要求产出实质性内容。
只改标签不写内容 = 空壳 S，会被质量自检系统检测出来。

### 问题 2：论文之间缺乏有意义的关联

关键词共现、引用图谱能告诉你"谁引用了谁"，但不能告诉你：

- 两篇论文是否共享同一个分子机制？
- 它们的结论是互相支持还是互相矛盾？
- A 的技术方法能否验证 B 的发现？

**LLS 的解决方案**：`gen_edges.py` v3.1 使用 5 个互补策略（显式引用、共享分子、
文本重叠、同病同法、隐藏轴），全部基于倒排索引实现 O(M) 复杂度，
生成有生物学意义的语义关联边。

### 问题 3：证据评估靠"感觉"

"这是 Nature 的论文，所以证据强"——期刊声望不等于证据质量。

**LLS 的解决方案**：内置标准化证据评分标准（`assets/data/evidence-rubric.json`），
基于偏倚风险、样本量、对照质量、可重复性等维度评分。
`enrich_paper_if()` 自动标注 IF 作为参考元数据，但不作为质量评分。

### 问题 4：知识无法积累和查询

读过的论文过几个月就忘了细节。要找"那篇关于 PD-L1 灰色地带的研究"时，
只能在 PDF 文件夹里靠文件名猜测。

**LLS 的解决方案**：所有分析结果以 NDJSON 格式持久化到知识图谱中。
`kg.py search` 支持全文检索，`build_network.py` 生成交互式网络图，
`export_citations.py` 导出 BibTeX——知识是可查询、可遍历、可复用的。

---

## 3. 核心能力

### 3.1 多源学术检索

一条命令搜索 PubMed、arXiv、bioRxiv、medRxiv、Crossref：

```bash
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 20
python scripts/literature_search.py arxiv "foundation model pathology" -y 2025
```

### 3.2 文献身份验证

在分析或引用之前，交叉验证文献身份：

```bash
python scripts/verify_citation.py --doi 10.xxxx/xxxxx
python scripts/verify_citation.py --pmid 12345678
```

检查撤稿、勘误、版本一致性。未验证的引用不会被加入知识图谱。

### 3.3 全文获取

按"最不脆弱"的顺序尝试获取全文：

1. **PubMed Central OA XML**（开放获取，无需认证）
2. **arXiv HTML**（预印本，curl 直接下载）
3. **bioRxiv/medRxiv Chrome CDP**（JavaScript 渲染页面，用户手动过验证）
4. **本地 PDF**（用户提供）

```bash
python scripts/fulltext_fetch.py --pmid 12345678
node scripts/extract_biorxiv_cdp.mjs --doi 10.1101/XXXX --port 9223
```

### 3.4 S 级 7 层深度分析

这是 LLS 最核心的能力（详见[§7](#7-核心方法论预览s-级-7-层解剖)）。每篇研究论文经过 7 层解剖：

| 层 | 名称 | 输出要求 |
|---|---|---|
| T1 | 文献基本档案 | 标题、作者、期刊（自动标注 IF）、研究设计、样本量 |
| T2 | 核心科学问题 | 1 个可证伪核心问题 + ≥5 个可检验子问题 |
| T3 | 主张-证据-综合链 | ≥5 条 CES 链，每条含具体数据（>20 字 evidence） |
| T4 | 分子机制级联 | ≥3 步因果链，精确到修饰位点，标注证据状态 |
| T5 | 隐藏组织轴 | ≥3 组观察-解读对，发现论文未明写的深层规律 |
| T6 | 概念创新图景 | 新概念、被推翻的旧观点、方法学突破、边界条件 |
| T7 | 跨文献关联 | ≥5 条有实质生物学描述（60-150 字）的跨论文关系 |

### 3.5 5 策略语义关联边生成

`gen_edges.py` v3.1 在论文间建立**有生物学意义**的关联：

| # | 策略 | 关系类型 | 数据基础 |
|---|---|---|---|
| 1 | 显式引用 | extends, accompanied_by | T7 cross_refs（LLM 深度分析结果） |
| 2 | 共享分子 | shares_molecules (≥2) | 90,125 基因符号，两段式精确匹配 |
| 2.5 | 文本重叠 | shares_topic (≥4) | 核心发现 bag-of-words 倒排索引 |
| 3 | 同病同法 | shares_disease_method | 疾病标签 × 方法标签交叉 |
| 4 | 隐藏轴 | shares_paradigm | T5 深层规律关键词共鸣 |
| 5 | 概念节点 | defines_concept | concepts.db → 论文 |

所有边都有自然语言描述说明 WHY。明确禁止 `same_journal`、`same_issue` 等无生物学意义的边。

### 3.6 自动影响因子标注

`kg_core.enrich_paper_if()` 从内置的 21,800 条 JCR 2024 期刊记录中自动匹配：

```
Nature → IF=50.5, Q1
Cell → IF=45.5, Q1
Science → IF=44.7, Q1
```

系统优先使用用户提供的期刊数据，再回退到打包数据。

### 3.7 10 维度质量自检

`selfcheck_knowledge_graph.py` 对知识图谱执行全面健康检查：

1. **文件清单**：全目录扫描，文件数/大小统计
2. **禁止残留**：Chrome CDP profile、cookie 文件、临时文件
3. **CDP 端口**：检查 9222/9223 是否仍打开
4. **DB 完整性**：NDJSON 解析错误、重复 ID、缺失必填字段
5. **S 级质量**：逐篇 7 项空壳检测（T2 子问题数、T3 链数、evidence 字数……）
6. **概念审计**：重复 ID、缺失 name
7. **边审计**：非法非生物关系、无描述边、孤儿边、自环、重复边
8. **全文缓存审计**：命名规范、过小文件、Cloudflare 残留
9. **边统计**：按策略的边数分布
10. **网络一致性**：交叉引用有效性验证

### 3.8 交互式可视化

```bash
python scripts/build_network.py
# → network.html（双击即开，纯离线）
```

- 节点颜色区分来源（PubMed/arXiv/bioRxiv/medRxiv/概念）
- 连线颜色区分关系类型
- 节点大小反映分析深度
- 支持力导向布局和交互探索

### 3.9 自动化监控

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json
```

监控器特性：可续（中断后继续）、批处理（避免内存溢出）、幂等（去重）、故障隔离。

---

## 4. 系统架构

### 4.1 数据流

```
研究问题
    │
    ▼
[1. 构架检索策略]   ← PICO/PECO 框架
    │
    ▼
[2. 多源检索]       ← PubMed / arXiv / bioRxiv / Crossref
    │
    ▼
[3. 身份验证]       ← 双源交叉验证 DOI/PMID
    │
    ▼
[4. 全文获取]       ← PMC XML / arXiv HTML / CDP 提取
    │
    ▼
[5. S 级 7 层分析]  ← LLM 深度推理（T1-T7）
    │
    ▼
[6. 证据评分]       ← 偏倚风险 / 样本量 / 可重复性
    │
    ▼
[7. 知识图谱持久化] ← NDJSON append
    │
    ├── [8a. 关联边生成]  gen_edges.py v3.1
    ├── [8b. 每日速报]    gen_digest.py
    ├── [8c. 网络可视化]  build_network.py
    └── [8d. 质量自检]    selfcheck_knowledge_graph.py
```

### 4.2 运行时工作区

`init_workspace.py` 创建的自包含工作区：

```
my-workspace/
├── workspace.json            ← 工作区元数据
├── papers.db                 ← NDJSON，每行一篇论文
├── concepts.db               ← NDJSON，每行一个概念节点
├── edges.db                  ← NDJSON，每行一条语义边
├── queries.db                ← NDJSON，检索记录
├── journal_metrics.db        ← 可选，用户自定义 IF 数据
├── data/                     ← 播种的基因/通路字典
│   ├── bioc_genes.json       ← 90,125 基因符号
│   └── kegg_pathways.json    ← 25,939 通路/GO 术语
├── fulltext/                 ← 下载的全文文档
├── fulltext_cache/           ← 全文缓存（gen_edges 使用）
├── daily_digest/             ← 每日速报 + 自检报告
├── reports/                  ← 生成的报告
├── exports/                  ← BibTeX 等导出
├── config/                   ← 监控配置
└── logs/                     ← 运行日志
```

所有 `.db` 文件为 **NDJSON 文本格式**（非 SQLite），可用任何文本编辑器打开，
`git diff` 友好，`head`/`grep`/`jq` 可直接查询。

---

## 5. 完整工作流

### 5.1 初始化

```bash
# 安装
git clone <repo-url> && cd literature-learning-suite
pip install -r scripts/requirements.txt

# 创建工作区
python scripts/init_workspace.py --root ./my-workspace

# 设置环境变量（脚本自动找到工作区）
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"
```

### 5.2 检索与入库

```bash
# 检索
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10 -o results.jsonl

# 标准化 + 验证
python scripts/normalize_records.py results.jsonl -o normalized.jsonl
python scripts/validate_records.py normalized.jsonl

# 入库（自动去重 + IF 标注）
python scripts/kg.py --root ./my-workspace add normalized.jsonl

# 查看统计
python scripts/kg.py --root ./my-workspace stats
```

### 5.3 深度分析后入库

通过 LLM 或人工完成 7 层分析后，将完整的 JSON 记录追加到 papers.db：

```bash
# 方式 1：通过 kg.py
python scripts/kg.py --root ./my-workspace add analyzed_papers.jsonl

# 方式 2：通过 Python API
python -c "
import sys; sys.path.insert(0, 'scripts')
from kg_core import add_paper, enrich_paper_if
paper = {'title': '...', 'pmid': '...', 'journal': 'Nature', 'tier2_core_question': '...', ...}
paper = enrich_paper_if(paper)
add_paper(paper)
"
```

### 5.4 关联边与可视化

```bash
# 关键：先清 .pyc 缓存
rm -rf scripts/__pycache__

# 生成语义关联边
python -B scripts/gen_edges.py

# 生成交互式网络图
python scripts/build_network.py

# 生成每日速报
python scripts/gen_digest.py
```

### 5.5 质量自检

```bash
python scripts/selfcheck_knowledge_graph.py
# 输出：daily_digest/selfcheck_YYYY-MM-DD.json + .md
```

### 5.6 导出

```bash
python scripts/export_citations.py ./my-workspace/papers.db --format bibtex -o library.bib
```

---

## 6. 安装与快速开始

### 环境要求

| 组件 | 最低版本 | 说明 |
|------|---------|------|
| Python | 3.10+ | 核心运行环境 |
| pip | 最新 | 包管理 |
| Node.js | v24+ | 仅 bioRxiv CDP 提取需要 |
| R + Bioconductor | — | 仅重新生成基因字典时需要 |

### 三步开始

```bash
# 1. 克隆并安装依赖
git clone <repo-url>
cd literature-learning-suite
pip install -r scripts/requirements.txt

# 2. 初始化工作区（播种基因/通路字典）
python scripts/init_workspace.py

# 3. 运行第一条检索
python scripts/literature_search.py pubmed "CRISPR screen cancer" -n 5
```

**无需额外配置**。打包的基因字典、通路术语和期刊指标开箱即用。
不需要安装数据库引擎、不需要申请 API 密钥（PubMed 直连）、不需要 Bioconductor。

---

## 7. 核心方法论预览：S 级 7 层解剖

> 完整方法论见 [GUIDE_ZH.md §4](GUIDE_ZH.md#4-核心方法论s级7层解剖协议)（中文，47KB）

### 为什么是 7 层

传统阅读：看标题 → 读摘要 → 标记 → 写两句。这种浅层处理**无法捕捉**：

- 论文的深层逻辑（实验为什么这样设计？隐含假设是什么？）
- 多条证据链的独立性（5 条 C-E-S 链要求从不同角度审视同一篇论文）
- 证据的真实强度（p 值 ≠ 效果量，统计显著 ≠ 临床显著）
- 跨论文的隐含关联（共享分子机制、隐藏范式共鸣）

7 层协议的每一层都要求产出**实质性内容**——只改标签不写内容 = 空壳 S，会被自检系统检出。

### 七层速览

| 层 | 名称 | 一句话 |
|---|---|---|
| T1 | 文献基本档案 | 这篇论文是什么？（可自动提取） |
| T2 | 核心科学问题 | 它要回答什么？拆成 ≥5 个可检验子问题 |
| T3 | 主张-证据-综合链 | 每一条主张有什么具体证据支撑？≥5 条独立链 |
| T4 | 分子机制级联 | 从触发到表型的完整因果链，精确到修饰位点 |
| T5 | 隐藏组织轴 | 论文**没有明写**的深层规律——你发现了什么？ |
| T6 | 概念创新图景 | 新概念？推翻的旧观点？方法学突破？边界在哪？ |
| T7 | 跨文献关联 | 和哪些论文有什么实质关系？为什么？≥5 条 |

### T3 示例

**好的 CES 链**：
> Claim：PD-L1 1-49% 亚组从新辅助化学免疫治疗中获益。
> Evidence：2847 例 NSCLC，SHAP 交互分析，EFS HR=0.52（95% CI 0.34-0.78），p=0.002，独立影像学评审。
> Synthesis：填补了 KEYNOTE-671 排除的 PD-L1 "灰色地带"人群的证据空白。
> Strength：★★★★（多中心 RCT，大样本，独立评审）
> Uncertain：检测平台一致性（22C3 vs 28-8），亚洲人群占比高。

**空壳 S（差）**：
> Claim：发现了新的生物标志物。
> Evidence：作者证明了。（≤20 字 → 空壳检测触发）

### T5 示例

> Observation：所有降维分析中，肿瘤边缘样本始终与核心样本分离为独立聚类。
> Interpretation：该研究**隐含地**将"边缘"而非"核心"定义为疾病决定性区室，
> 这解释了为什么核心基因特征的预后价值反而更低。

T5 考察的是你**自己**的发现能力——不是复述 Discussion。

### 空壳检测

写后立即运行验证：

```bash
python -c "
import json
with open('./my-workspace/papers.db', encoding='utf-8-sig') as f:
    for line in f:
        p = json.loads(line)
        if p.get('analysis_tier') != 'S': continue
        chains = p.get('tier3_ces_chains', [])
        empty = sum(1 for c in chains if len(c.get('evidence','')) <= 20)
        print(f'{p[\"id\"]}: T2={len(p.get(\"tier2_subquestions\",[]))} '
              f'T3={len(chains)} empty_ev={empty} '
              f'T5={len(p.get(\"tier5_hidden_axis\",[]))} '
              f'T7={len(p.get(\"tier7_cross_refs\",[]))}')
"
```

---

## 8. 工具链概览

### 检索与发现
| 脚本 | 功能 |
|------|------|
| `literature_search.py` | 多源检索（PubMed/arXiv/Crossref/bioRxiv） |
| `search_arxiv.py` | arXiv 专项（按作者/分类/日期过滤） |
| `download_biorxiv_api.py` | bioRxiv/medRxiv 批量下载 |

### 验证与标准化
| 脚本 | 功能 |
|------|------|
| `verify_citation.py` | DOI/PMID/arXiv 交叉验证 |
| `normalize_records.py` | 标准化 + 去重（PMID > DOI > 标题+年份） |
| `validate_records.py` | JSON Schema 验证 |

### 全文获取
| 脚本 | 功能 |
|------|------|
| `fulltext_fetch.py` | 统一下载（PMC XML + arXiv HTML） |
| `extract_pymupdf.py` | PDF 文本/表格（PyMuPDF） |
| `extract_marker.py` | PDF OCR（Marker） |
| `extract_biorxiv_cdp.mjs` | Chrome CDP 全文提取（Node.js） |

### 知识图谱
| 脚本 | 功能 |
|------|------|
| `kg.py` | CLI：add / stats / search / audit |
| `kg_core.py` | 核心库：CRUD + 自动 IF 标注 |

### 分析与输出
| 脚本 | 功能 |
|------|------|
| `gen_edges.py` | 语义关联边生成器 v3.1（5 策略） |
| `gen_digest.py` | 每日速报 |
| `build_network.py` | 交互式力导向 HTML 网络图 |
| `selfcheck_knowledge_graph.py` | 10 维度质量自检 |
| `export_citations.py` | BibTeX / CSL 导出 |

### 自动化
| 脚本 | 功能 |
|------|------|
| `monitor.py` | 可续批监控器 |
| `init_workspace.py` | 工作区初始化 + 数据播种 |

---

## 9. MCP 集成（可选增强）

LLS 在无 MCP 的情况下完全可用，但接入 MCP 后 Agent 可以直接调用检索工具，
无需手动运行 CLI 脚本，工作流更流畅。

### 支持的 MCP 服务器

| MCP 服务器 | 用途 | 安装 |
|-----------|------|------|
| **PubMed MCP** | 检索 PubMed + 获取 PMC 全文 | `pip install mcp-simple-pubmed` |
| **arXiv MCP** | 搜索/下载/引用图谱/主题监控 | `pip install arxiv-mcp-server` |
| **Fetch MCP** | 抓取公开网页 | `pip install mcp-server-fetch` |
| **Playwright MCP** | JS 渲染页面（bioRxiv 正文） | `npx @playwright/mcp@latest` |

### 快速配置

将 `assets/templates/mcp-servers.yaml` 的内容合并到你的 Agent 宿主配置中：

```yaml
mcp_servers:
  pubmed:
    command: mcp-simple-pubmed
    env:
      PUBMED_EMAIL: "your@email.edu"     # 必填
  arxiv:
    command: arxiv-mcp-server
    args:
      - --storage-path
      - "./arxiv-storage"
  fetch:
    command: mcp-server-fetch
```

### 工具回退

| 需求 | 优先（MCP） | 回退（CLI） |
|------|------------|------------|
| PubMed 检索 | `search_pubmed()` | `literature_search.py pubmed` |
| arXiv 检索 | `search_papers()` | `literature_search.py arxiv` |
| PMC 全文 | `get_paper_fulltext()` | `fulltext_fetch.py --pmid` |
| JS 渲染页 | Playwright MCP | `extract_biorxiv_cdp.mjs` |

📖 完整设置指南：[references/mcp-setup.md](references/mcp-setup.md)
📖 工具路由策略：[references/mcp-and-tool-routing.md](references/mcp-and-tool-routing.md)

---

## 10. 文档导航

### 使用指南（7 种语言）

| 语言 | 文件 | 说明 |
|------|------|------|
| 🇨🇳 **中文** | [**GUIDE_ZH.md**](GUIDE_ZH.md) | **主语言**，47KB，含完整方法论、全部示例、陷阱库、FAQ |
| 🇬🇧 English | [GUIDE_EN.md](GUIDE_EN.md) | 22KB，完整英文指南 |
| 🇩🇪 Deutsch | [GUIDE_DE.md](GUIDE_DE.md) | 54KB，完整德文指南 |
| 🇯🇵 日本語 | [GUIDE_JA.md](GUIDE_JA.md) | 58KB，完整日文指南 |
| 🇰🇷 한국어 | [GUIDE_KO.md](GUIDE_KO.md) | 55KB，完整韩文指南 |
| 🇪🇸 Español | [GUIDE_ES.md](GUIDE_ES.md) | 28KB，完整西班牙文指南 |
| 🇸🇦 العربية | [GUIDE_AR.md](GUIDE_AR.md) | 32KB，完整阿拉伯文指南 |

### 操作协议与参考文档

| 文件 | 说明 |
|------|------|
| [SKILL.md](SKILL.md) | AI Agent 操作协议——Agent 加载此文件了解完整操作流程 |
| [references/](references/) | 25+ 篇方法论文档（数据模型、深层分析、边生成、全文获取、排障等） |
| [assets/schemas/](assets/schemas/) | JSON Schema（paper/concept/edge/query/monitor/workspace/journal） |
| [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md) | 第三方数据归属说明 |

---

## 10. 打包数据

| 数据文件 | 内容 | 记录数 | 来源 |
|---------|------|--------|------|
| `assets/data/bioc_genes.json` | 人 + 小鼠基因符号 | 90,125 | `org.Hs.eg.db` + `org.Mm.eg.db` (Bioconductor) |
| `assets/data/kegg_pathways.json` | KEGG 通路 + GO BP 术语 | 25,939 | KEGGREST + `GO.db` |
| `assets/data/journal_metrics_2024.json` | JCR 2024 期刊 IF、分区、分类 | 21,800 | Clarivate Journal Citation Reports |
| `assets/data/evidence-rubric.json` | 证据评分标准 | — | EBM 方法论 |
| `assets/data/relation-ontology.json` | 关联边关系类型本体 | — | 文献 KG 设计 |
| `assets/data/study-designs.json` | 研究设计分类 | — | Cochrane / CEBM |
| `assets/data/search-query-packs.json` | 预配置检索模板 | — | 人工整理 |

所有数据附带 **SHA-256 校验**（`assets/data/data-manifest.json`），可通过 `python scripts/check_assets.py` 验证。

基因本体 (GO) 内容遵循 GO 归属和许可条款。KEGG 内容的再分发条款因使用场景而异。JCR 期刊指标来自 Clarivate Journal Citation Reports 2024 版。再分发修改版或重新生成的副本前，请查阅相应的上游条款。

---

## 11. 许可证

### CC BY-NC-SA 4.0 — 署名-非商业性使用-相同方式共享

**你可以自由地：**
- **分享** — 在任何媒介以任何形式复制、分发本作品
- **改编** — 修改、转换或以本作品为基础进行创作

**但必须遵守以下条件：**
- **署名** — 必须给出适当的署名，提供指向本许可证的链接，同时标明是否对原始作品作了修改
- **非商业性使用** — 不得将本作品用于**商业目的**
- **相同方式共享** — 如果你再混合、转换或基于本作品进行创作，**必须基于与原先相同的许可证**分发你的贡献作品

**这意味着：**

| 场景 | 是否允许 |
|------|---------|
| 个人学术研究使用 | ✅ 允许 |
| 大学/研究所内部使用 | ✅ 允许 |
| 修改后用于自己的论文分析 | ✅ 允许（需署名） |
| 将修改版发布到 GitHub | ✅ 允许（必须以 CC BY-NC-SA 4.0 开源） |
| 用于公司内部研究 | ✅ 允许（非商业） |
| 打包为 SaaS 产品的一部分出售 | ❌ 禁止 |
| 嵌入商业软件分发 | ❌ 禁止 |
| 修改后闭源分发 | ❌ 禁止（违反 ShareAlike） |
| 用于商业培训课程 | ❌ 禁止 |

**二创产品必须开源**：如果你基于 LLS 创建了衍生作品（例如修改版、扩展版、包装版），你必须以 **CC BY-NC-SA 4.0 或兼容许可证** 开源你的整个衍生作品。

全文见 [LICENSE](LICENSE)（含完整法律文本和人话摘要）。

---

## 12. 目录结构

```
literature-learning-suite/
│
├── GUIDE_ZH.md                    ← 中文完整指南（主语言，47KB）
├── GUIDE_EN.md                    ← English guide (22KB)
├── GUIDE_DE.md                    ← Deutsche Anleitung (54KB)
├── GUIDE_JA.md                    ← 日本語ガイド (58KB)
├── GUIDE_KO.md                    ← 한국어 가이드 (55KB)
├── SKILL.md                       ← AI Agent 操作协议
├── README.md                      ← 本文件
├── LICENSE                        ← CC BY-NC-SA 4.0 完整文本
├── THIRD_PARTY_DATA.md            ← 第三方数据归属
├── VERSION                        ← 语义版本号
│
├── scripts/                       ← 核心工具链
│   ├── init_workspace.py          ← 工作区初始化
│   ├── literature_search.py       ← 多源检索
│   ├── search_arxiv.py            ← arXiv 专项检索
│   ├── download_biorxiv_api.py    ← bioRxiv API 下载
│   ├── verify_citation.py         ← 文献验证
│   ├── normalize_records.py       ← 标准化去重
│   ├── validate_records.py        ← JSON Schema 验证
│   ├── fulltext_fetch.py          ← 全文统一下载
│   ├── extract_pymupdf.py         ← PDF 提取 (PyMuPDF)
│   ├── extract_marker.py          ← PDF OCR (Marker)
│   ├── extract_biorxiv_cdp.mjs    ← Chrome CDP 全文 (Node.js)
│   ├── biorxiv_chrome_cdp_launcher.bat  ← CDP 启动器 (Win)
│   ├── kg.py                      ← KG CLI (add/stats/search/audit)
│   ├── kg_core.py                 ← KG 核心库 + IF 标注
│   ├── ll_common.py               ← 共享工具
│   ├── workspace_paths.py         ← 路径解析
│   ├── gen_edges.py               ← 关联边生成 v3.1
│   ├── gen_digest.py              ← 每日速报
│   ├── build_network.py           ← 交互式网络图
│   ├── selfcheck_knowledge_graph.py  ← 质量自检
│   ├── export_citations.py        ← BibTeX 导出
│   ├── export_bioc_genes.R        ← 基因字典再生 (R)
│   ├── journal_metrics.py         ← 期刊指标导入
│   ├── monitor.py                 ← 批监控器
│   ├── check_assets.py            ← 资产校验
│   └── requirements.txt           ← Python 依赖
│
├── assets/
│   ├── data/                      ← 打包参考数据（7 文件，附 SHA-256）
│   ├── schemas/                   ← JSON Schema（7 种记录类型）
│   └── templates/                 ← 记录模板 + 配置模板
│
├── references/                    ← 方法论与操作协议（25+ 篇）
│
└── tests/                         ← 单元测试（18 个，全部通过）
    └── fixtures/                  ← 合成测试数据
```

---

> **Literature Learning Suite v1.3.0** · CC BY-NC-SA 4.0
>
> 从研究问题到知识图谱。不是管理文献，是理解文献。
