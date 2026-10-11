# Literature Learning Suite

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.3.0-green.svg)](VERSION)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](scripts/requirements.txt)

学术文献的深度分析与知识图谱构建工具，供人和 AI Agent 使用。检索、全文获取、7 层结构化解剖、语义关联边生成和每日速报在同一套工作区里完成，分析结果以 NDJSON 知识图谱保存在本地，可以全文检索、图遍历和导出。

## 它解决什么问题

- 阅读停留在表面：只看标题和摘要，抓不到论文的深层逻辑、隐含假设，以及统计显著和临床显著的区别。7 层协议为每层写死了产出要求，只改标签的条目会被判为空壳。
- 论文之间缺少有意义的关联：关键词共现和引用图谱只说明谁引用了谁。5 种策略基于共享分子、文本重叠、同病同法和隐藏轴生成有生物学含义的语义边。
- 证据评估靠感觉：期刊声望不等于证据质量。内置评分标准从偏倚风险、样本量、对照质量和可重复性打分，影响因子只作参考元数据，不计入评分。
- 知识无法积累：读过的论文过几个月只剩文件名。分析结果持久化为 NDJSON，支持全文检索、图遍历和 BibTeX 导出。

## 安装

| 组件 | 最低版本 | 说明 |
|------|---------|------|
| Python | 3.10+ | 核心运行环境 |
| pip | 最新 | 包管理 |
| Node.js | v24+ | 仅 bioRxiv CDP 提取需要 |
| R + Bioconductor | | 仅重新生成基因字典时需要 |

```bash
# 1. 克隆并安装依赖
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/literature-learning-suite
pip install -r scripts/requirements.txt

# 2. 初始化工作区（播种基因与通路字典）
python scripts/init_workspace.py

# 3. 运行第一条检索
python scripts/literature_search.py pubmed "CRISPR screen cancer" -n 5
```

装进 Hermes Agent 后可以直接用对话驱动：

```bash
cp -r OpenBioMed-Lab/literature-learning-suite ~/.hermes/skills/
```

`init_workspace.py` 创建自包含工作区，基因字典、通路术语和期刊指标随包提供，不需要安装数据库引擎、申请 API 密钥（PubMed 直连）或使用 Bioconductor。脚本通过 `LITERATURE_KG_ROOT` 找工作区：

```bash
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"
```

## 用法

### 检索与入库

```bash
# 多源检索
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10 -o results.jsonl
python scripts/literature_search.py arxiv "foundation model pathology" -y 2025

# 交叉验证文献身份（检查撤稿、勘误、版本一致性）
python scripts/verify_citation.py --doi 10.xxxx/xxxxx
python scripts/verify_citation.py --pmid 12345678

# 标准化 + 校验 + 入库（自动去重并标注影响因子）
python scripts/normalize_records.py results.jsonl -o normalized.jsonl
python scripts/validate_records.py normalized.jsonl
python scripts/kg.py --root ./my-workspace add normalized.jsonl
python scripts/kg.py --root ./my-workspace stats
```

未通过身份验证的引用不会被加入知识图谱。

### 全文获取

按最不脆弱的顺序尝试，前面成功就不再往下走：

1. PubMed Central OA XML（开放获取，无需认证）
2. arXiv HTML（预印本，直接下载）
3. bioRxiv/medRxiv 的 Chrome CDP（JavaScript 渲染页面，用户手动过验证）
4. 本地 PDF（用户提供）

```bash
python scripts/fulltext_fetch.py --pmid 12345678
node scripts/extract_biorxiv_cdp.mjs --doi 10.1101/XXXX --port 9223
```

### 深度分析

经过 7 层解剖的论文用 `kg.py` 或 Python API 写入知识图谱：

```bash
python scripts/kg.py --root ./my-workspace add analyzed_papers.jsonl
```

```python
import sys; sys.path.insert(0, 'scripts')
from kg_core import add_paper, enrich_paper_if

paper = {'title': '...', 'pmid': '...', 'journal': 'Nature', 'tier2_core_question': '...'}
paper = enrich_paper_if(paper)
add_paper(paper)
```

### 关联边与可视化

```bash
# 生成语义关联边（先清 .pyc 缓存）
rm -rf scripts/__pycache__
python -B scripts/gen_edges.py

# 交互式网络图与每日速报
python scripts/build_network.py
python scripts/gen_digest.py
```

`build_network.py` 生成单文件 `network.html`，双击即开，离线可用。节点颜色区分来源（PubMed/arXiv/bioRxiv/medRxiv/概念），连线颜色区分关系类型，节点大小反映分析深度，支持力导向布局和交互探索。

### 质量自检

```bash
python scripts/selfcheck_knowledge_graph.py
# 输出：daily_digest/selfcheck_YYYY-MM-DD.json + .md
```

### 导出

```bash
python scripts/export_citations.py ./my-workspace/papers.db --format bibtex -o library.bib
```

### 自动化监控

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json
```

监控器可续跑（中断后继续）、批次处理（避免内存溢出）、幂等（去重）、故障隔离。

## 核心方法：7 层结构化解剖

一篇论文经过 7 层解剖，每层都要求产出实质内容。只改标签不写内容会被质量自检判为空壳。协议全文与判分细则见 [GUIDE_ZH.md](GUIDE_ZH.md#4-核心方法论s级7层解剖协议)（中文）。

| 层 | 名称 | 输出要求 |
|---|---|---|
| T1 | 文献基本档案 | 标题、作者、期刊（自动标注影响因子）、研究设计、样本量 |
| T2 | 核心科学问题 | 1 个可证伪核心问题加至少 5 个可检验子问题 |
| T3 | 主张-证据-综合链 | 至少 5 条 CES 链，每条含具体数据（evidence 超过 20 字） |
| T4 | 分子机制级联 | 至少 3 步因果链，精确到修饰位点，标注证据状态 |
| T5 | 隐藏组织轴 | 至少 3 组观察-解读对，找出论文没有明写的深层规律 |
| T6 | 概念创新图景 | 新概念、被推翻的旧观点、方法学突破、边界条件 |
| T7 | 跨文献关联 | 至少 5 条有实质生物学描述的跨论文关系（60 到 150 字） |

传统阅读只做标题、摘要、标记和两句总结，抓不到论文的深层逻辑、多条证据链的独立性，以及统计显著和临床显著的区别。7 层协议把这些拆成可检查的产出。

T3 的一条合格 CES 链：

> Claim：PD-L1 1-49% 亚组从新辅助化学免疫治疗中获益。
> Evidence：2847 例 NSCLC，SHAP 交互分析，EFS HR=0.52（95% CI 0.34-0.78），p=0.002，独立影像学评审。
> Synthesis：填补了 KEYNOTE-671 排除的 PD-L1 "灰色地带"人群的证据空白。
> Strength：多中心 RCT，大样本，独立评审。
> Uncertain：检测平台一致性（22C3 与 28-8），亚洲人群占比高。

不合格的例子（evidence 不超过 20 字即触发空壳检测）：

> Claim：发现了新的生物标志物。
> Evidence：作者证明了。

T5 考察自己的发现能力，而不是复述 Discussion：

> Observation：所有降维分析中，肿瘤边缘样本始终与核心样本分离为独立聚类。
> Interpretation：该研究隐含地把"边缘"而不是"核心"定义为疾病决定性区室，这解释了为什么核心基因特征的预后价值反而更低。

写完立即验证：

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

### 语义关联边

`gen_edges.py` v3.1 用 5 个策略建立论文之间的关联，全部基于倒排索引，复杂度 O(M)：

| # | 策略 | 关系类型 | 数据基础 |
|---|---|---|---|
| 1 | 显式引用 | extends, accompanied_by | T7 cross_refs（LLM 深度分析结果） |
| 2 | 共享分子 | shares_molecules（至少 2 个） | 90,125 个基因符号，两段式精确匹配 |
| 2.5 | 文本重叠 | shares_topic（至少 4 个） | 核心发现 bag-of-words 倒排索引 |
| 3 | 同病同法 | shares_disease_method | 疾病标签与方向标签交叉 |
| 4 | 隐藏轴 | shares_paradigm | T5 深层规律关键词共鸣 |
| 5 | 概念节点 | defines_concept | concepts.db 到论文 |

每条边都有自然语言描述说明为什么建立。`same_journal`、`same_issue` 这类没有生物学含义的边被明确禁止。

### 影响因子标注

`kg_core.enrich_paper_if()` 从随包的 21,800 条 JCR 2024 期刊记录中匹配：

```text
Nature  -> IF=50.5, Q1
Cell    -> IF=45.5, Q1
Science -> IF=44.7, Q1
```

用户提供期刊数据时优先使用用户数据，再回退到随包数据。影响因子只作参考元数据，不计入证据评分。

### 10 维度质量自检

`selfcheck_knowledge_graph.py` 检查：

1. 文件清单：全目录扫描，统计文件数和大小
2. 禁止残留：Chrome CDP profile、cookie 文件、临时文件
3. CDP 端口：检查 9222 与 9223 是否仍然打开
4. DB 完整性：NDJSON 解析错误、重复 ID、缺失必填字段
5. S 级质量：逐篇做 7 项空壳检测（T2 子问题数、T3 链数、evidence 字数等）
6. 概念审计：重复 ID、缺失 name
7. 边审计：非法非生物关系、无描述边、孤儿边、自环、重复边
8. 全文缓存审计：命名规范、过小文件、Cloudflare 残留
9. 边统计：按策略统计边数分布
10. 网络一致性：交叉引用有效性验证

### 数据流

1. 构架检索策略（PICO/PECO 框架）
2. 多源检索（PubMed、arXiv、bioRxiv、Crossref）
3. 身份验证（DOI/PMID 双源交叉验证）
4. 全文获取（PMC XML、arXiv HTML、CDP 提取）
5. S 级 7 层分析（LLM 深度推理，T1 到 T7）
6. 证据评分（偏倚风险、样本量、可重复性）
7. 知识图谱持久化（NDJSON 追加）
8. 关联边生成、每日速报、网络可视化、质量自检

### 工作区结构

`init_workspace.py` 创建的自包含工作区：

```text
my-workspace/
├── workspace.json            工作区元数据
├── papers.db                 NDJSON，每行一篇论文
├── concepts.db               NDJSON，每行一个概念节点
├── edges.db                  NDJSON，每行一条语义边
├── queries.db                NDJSON，检索记录
├── journal_metrics.db        可选，用户自定义影响因子数据
├── data/                     播种的基因与通路字典
│   ├── bioc_genes.json       90,125 个基因符号
│   └── kegg_pathways.json    25,939 条通路与 GO 术语
├── fulltext/                 下载的全文文档
├── fulltext_cache/           全文缓存（gen_edges 使用）
├── daily_digest/             每日速报与自检报告
├── reports/                  生成的报告
├── exports/                  BibTeX 等导出
├── config/                   监控配置
└── logs/                     运行日志
```

所有 `.db` 文件是 NDJSON 文本，不是 SQLite，可以用任何文本编辑器打开，`git diff` 友好，`head`、`grep`、`jq` 可以直接查询。

## 命令与脚本

检索与发现：

| 脚本 | 功能 |
|------|------|
| `literature_search.py` | 多源检索（PubMed、arXiv、Crossref、bioRxiv） |
| `search_arxiv.py` | arXiv 专项检索（按作者、分类、日期过滤） |
| `download_biorxiv_api.py` | bioRxiv 与 medRxiv 批量下载 |

验证与标准化：

| 脚本 | 功能 |
|------|------|
| `verify_citation.py` | DOI、PMID、arXiv 交叉验证 |
| `normalize_records.py` | 标准化与去重（PMID 优先于 DOI，再退到标题加年份） |
| `validate_records.py` | JSON Schema 验证 |

全文获取：

| 脚本 | 功能 |
|------|------|
| `fulltext_fetch.py` | 统一下载（PMC XML 与 arXiv HTML） |
| `extract_pymupdf.py` | PDF 文本与表格（PyMuPDF） |
| `extract_marker.py` | PDF OCR（Marker） |
| `extract_biorxiv_cdp.mjs` | Chrome CDP 全文提取（Node.js） |
| `biorxiv_chrome_cdp_launcher.bat` | CDP 启动器（Windows） |

知识图谱：

| 脚本 | 功能 |
|------|------|
| `kg.py` | 命令行入口：add、stats、search（全文检索）、audit |
| `kg_core.py` | 核心库：增删改查与影响因子标注 |
| `ll_common.py`、`workspace_paths.py` | 共享工具与路径解析 |

分析与输出：

| 脚本 | 功能 |
|------|------|
| `gen_edges.py` | 语义关联边生成 v3.1（5 策略） |
| `gen_digest.py` | 每日速报 |
| `build_network.py` | 交互式力导向网络图 |
| `selfcheck_knowledge_graph.py` | 10 维度质量自检 |
| `export_citations.py` | BibTeX 与 CSL 导出 |

维护：

| 脚本 | 功能 |
|------|------|
| `monitor.py` | 可续跑批监控器 |
| `init_workspace.py` | 工作区初始化与数据播种 |
| `check_assets.py` | 随包数据校验（SHA-256） |
| `journal_metrics.py` | 期刊指标导入 |
| `export_bioc_genes.R` | 基因字典再生（R 与 Bioconductor） |
| `requirements.txt` | Python 依赖 |

## 打包数据

| 数据文件 | 内容 | 记录数 | 来源 |
|---------|------|--------|------|
| `assets/data/bioc_genes.json` | 人与小鼠基因符号 | 90,125 | `org.Hs.eg.db`、`org.Mm.eg.db`（Bioconductor） |
| `assets/data/kegg_pathways.json` | KEGG 通路与 GO BP 术语 | 25,939 | KEGGREST、`GO.db` |
| `assets/data/journal_metrics_2024.json` | JCR 2024 期刊影响因子、分区、分类 | 21,800 | Clarivate Journal Citation Reports |
| `assets/data/evidence-rubric.json` | 证据评分标准 | | EBM 方法论 |
| `assets/data/relation-ontology.json` | 关联边关系类型本体 | | 文献知识图谱设计 |
| `assets/data/study-designs.json` | 研究设计分类 | | Cochrane、CEBM |
| `assets/data/search-query-packs.json` | 预配置检索模板 | | 人工整理 |

所有数据附带 SHA-256 校验（`assets/data/data-manifest.json`），可用 `python scripts/check_assets.py` 验证。

基因本体（GO）内容遵循 GO 的归属和许可条款。KEGG 内容的再分发条款因使用场景而异。JCR 期刊指标来自 Clarivate Journal Citation Reports 2024 版。再分发修改版或重新生成的副本前，请查阅相应的上游条款。

## MCP 集成（可选）

本项目在没有 MCP 的情况下完全可用。接入后 Agent 可以直接调用检索工具，不必手动运行命令行脚本。

| MCP 服务器 | 用途 | 安装 |
|-----------|------|------|
| PubMed MCP | 检索 PubMed 并获取 PMC 全文 | `pip install mcp-simple-pubmed` |
| arXiv MCP | 搜索、下载、引用图谱、主题监控 | `pip install arxiv-mcp-server` |
| Fetch MCP | 抓取公开网页 | `pip install mcp-server-fetch` |
| Playwright MCP | JavaScript 渲染页面（bioRxiv 正文） | `npx @playwright/mcp@latest` |

把 `assets/templates/mcp-servers.yaml` 的内容合并到 Agent 宿主配置：

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

工具回退：

| 需求 | 优先（MCP） | 回退（命令行） |
|------|------------|------------|
| PubMed 检索 | `search_pubmed()` | `literature_search.py pubmed` |
| arXiv 检索 | `search_papers()` | `literature_search.py arxiv` |
| PMC 全文 | `get_paper_fulltext()` | `fulltext_fetch.py --pmid` |
| JavaScript 渲染页 | Playwright MCP | `extract_biorxiv_cdp.mjs` |

完整设置指南见 [references/mcp-setup.md](references/mcp-setup.md)，工具路由策略见 [references/mcp-and-tool-routing.md](references/mcp-and-tool-routing.md)。

## 文档

| 语言 | 文件 | 说明 |
|------|------|------|
| 中文 | [GUIDE_ZH.md](GUIDE_ZH.md) | 主语言，47KB，含完整方法论、全部示例、陷阱库和 FAQ |
| English | [GUIDE_EN.md](GUIDE_EN.md) | 22KB |
| Deutsch | [GUIDE_DE.md](GUIDE_DE.md) | 54KB |
| 日本語 | [GUIDE_JA.md](GUIDE_JA.md) | 58KB |
| 한국어 | [GUIDE_KO.md](GUIDE_KO.md) | 55KB |
| Español | [GUIDE_ES.md](GUIDE_ES.md) | 28KB |
| العربية | [GUIDE_AR.md](GUIDE_AR.md) | 32KB |

| 文件 | 说明 |
|------|------|
| [SKILL.md](SKILL.md) | AI Agent 操作协议，Agent 加载此文件了解完整操作流程 |
| [references/](references/) | 25 篇以上方法论文档（数据模型、深层分析、边生成、全文获取、排障等） |
| [assets/schemas/](assets/schemas/) | JSON Schema（paper、concept、edge、query、monitor、workspace、journal） |
| [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md) | 第三方数据归属说明 |

## 许可证

CC BY-NC-SA 4.0，署名-非商业性使用-相同方式共享。完整法律文本和人话摘要见 [LICENSE](LICENSE)。

你可以自由地分享（在任何媒介以任何形式复制、分发本作品）和改编（修改、转换或以本作品为基础进行创作），但必须遵守以下条件：

- 署名：给出适当的署名，提供指向本许可证的链接，并标明是否作了修改
- 非商业性使用：不得将本作品用于商业目的
- 相同方式共享：再混合、转换或基于本作品创作时，必须以相同的许可证分发你的贡献作品

| 场景 | 是否允许 |
|------|---------|
| 个人学术研究使用 | 允许 |
| 大学或研究所内部使用 | 允许 |
| 修改后用于自己的论文分析 | 允许（需署名） |
| 将修改版发布到 GitHub | 允许（必须以 CC BY-NC-SA 4.0 开源） |
| 用于公司内部研究 | 允许（非商业） |
| 打包为 SaaS 产品的一部分出售 | 禁止 |
| 嵌入商业软件分发 | 禁止 |
| 修改后闭源分发 | 禁止（违反相同方式共享条款） |
| 用于商业培训课程 | 禁止 |

基于本项目创建衍生作品（修改版、扩展版、包装版）时，必须以 CC BY-NC-SA 4.0 或兼容许可证开源整个衍生作品。

## 署名

OpenBioMed Lab

| 姓名 | 机构 |
|------|------|
| Jinghua Gu | 安徽医科大学 · [2040519464@qq.com](mailto:2040519464@qq.com) · [ORCID: 0009-0000-8691-1312](https://orcid.org/0009-0000-8691-1312) |
| Shuyan Sheng | 安徽医科大学 |
| Huake Cao | 安徽医科大学 |
| Conghan Li | 安徽医科大学 |
| Taiyu Shi | 安徽医科大学 |

## 目录结构

```text
literature-learning-suite/
├── GUIDE_ZH.md                   中文完整指南（主语言，47KB）
├── GUIDE_EN.md                   English guide (22KB)
├── GUIDE_DE.md                   Deutsche Anleitung (54KB)
├── GUIDE_JA.md                   日本語ガイド (58KB)
├── GUIDE_KO.md                   한국어 가이드 (55KB)
├── GUIDE_ES.md                   Guía en español (28KB)
├── GUIDE_AR.md                   الدليل العربي (32KB)
├── SKILL.md                      AI Agent 操作协议
├── README.md                     本文件
├── LICENSE                       CC BY-NC-SA 4.0 完整文本
├── THIRD_PARTY_DATA.md           第三方数据归属
├── VERSION                       语义版本号
├── scripts/                      工具链，见「命令与脚本」
├── assets/
│   ├── data/                     随包参考数据（7 个文件，附 SHA-256）
│   ├── schemas/                  JSON Schema（7 种记录类型）
│   └── templates/                记录模板与配置模板
└── references/                   方法论与操作协议（25 篇以上）
```
