# 学术 API 连通性实测记录

> 环境：Windows 10 Pro @ 中国，Clash 代理 7890，节点 美国西雅图
> 测试日期：2026-06-08（更新：bioRxiv JATS/PDF 403 + JS渲染 + PubMed需unset代理）

## PubMed (eutils.ncbi.nlm.nih.gov)

**连通性**: ✅ 直连可用
**响应时间**: ~780ms
**推荐**: 首选检索源，**不需要代理**。走代理反而超时/连接失败——需 `unset http_proxy https_proxy`

### ESearch — 搜索
```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer+immunotherapy&retmax=3"
```

### EFetch — 获取全文摘要（批量，上限~200 IDs）
```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=PID1,PID2&retmode=xml&rettype=abstract"
```

### 日期过滤
```bash
BEFORE7D=$(date -d "7 days ago" +%Y/%m/%d)
TODAY=$(date +%Y/%m/%d)
term=...AND+(${BEFORE7D}[EDAT]:${TODAY}[EDAT])
```

## arXiv (export.arxiv.org)

**连通性**: ✅ 通过 Clash 代理可用
**关键发现**: 
- API 端点是 **HTTP**（不是 HTTPS）
- 请求 HTTP → 301 重定向到 HTTPS → **必须使用 `curl -L`** 跟随重定向
- **必须通过代理** `export http_proxy=http://127.0.0.1:7890`

```bash
export http_proxy=http://127.0.0.1:7890
curl -sL "http://export.arxiv.org/api/query?search_query=cat:q-bio.GN&sortBy=submittedDate&sortOrder=descending&max_results=20"
```

### 生物学相关类别
| 类别 | 字段 |
|------|------|
| `q-bio.GN` | Genomics |
| `q-bio.QM` | Quantitative Methods |
| `q-bio.BM` | Biomolecules |
| `q-bio.CB` | Cell Behavior |
| `q-bio.TO` | Tissues and Organs |
| `q-bio.MN` | Molecular Networks |
| `cs.CV` | Computer Vision (病理/医学图像) |
| `cs.LG` | Machine Learning |
| `cs.AI` | Artificial Intelligence |
| `cs.CL` | Computation and Language (生物NLP) |

## bioRxiv API

**连通性**: ✅ API 可用（通过代理获取元数据/摘要）
**全文访问**: ❌ JATS XML (`.source.xml`) 和 PDF (`.full.pdf`) 端点直接返回 403（2026-06-08 实测）。HTML 页面正文是 JS 动态渲染（React SPA），curl 即使携带有效 `cf_clearance` cookie 也拿不到实质内容
**返回**: JSON `{collection: [{doi, title, authors, abstract, category, date, ...}], messages: [{total: N}]}`
**策略**: API 获取 300-500 词摘要 → LLM 推理分析 → S 级 7 层解剖

```bash
export http_proxy=http://127.0.0.1:7890
curl -s "https://api.biorxiv.org/details/biorxiv/2026-06-07/2026-06-07/0/50"
```

## medRxiv API

**连通性**: ✅ API 可用（通过代理获取元数据/摘要）
**全文访问**: ❌ 同 bioRxiv——JATS XML 403、PDF 403、HTML JS 渲染
**结构**: 与 bioRxiv 相同
**策略**: API 获取摘要 → LLM 推理

```bash
export http_proxy=http://127.0.0.1:7890
curl -s "https://api.medrxiv.org/details/medrxiv/2026-06-07/2026-06-07/0/50"
```

### bioRxiv/medRxiv 全文获取唯一路径

**Playwright MCP**（真实 Chrome 浏览器 → 通过 CF → 执行 JS → 提取正文）。
详见 `references/preprint-fulltext.md`。

或在本地 Chrome 手动打开论文页面 → F12 → 复制正文 → 保存到 `fulltext_cache/`。

## PubMed Central (PMC) 全文

**连通性**: ⚠️ 仅 OA 期刊可用
- Cell/Nature/Science/Lancet/JAMA/Nature Medicine: **付费墙，PMC 无全文**
- OA 期刊（eLife/PLOS/BMC/Frontiers/Scientific Reports/Nature Communications）: PMC 可能有全文
**工具**: `mcp_pubmed_get_paper_fulltext(pmid="...")`

## 全文访问总结（2026-06-08 更新）

| 来源 | 元数据 | 全文 | 工具 | 备注 |
|------|--------|------|------|------|
| PubMed | ✅ 直连 | ⚠️ 仅OA | ESearch+EFetch | 国内直连无代理，走代理反而超时 |
| bioRxiv API | ✅ 代理 | — | JSON 摘要 | 300-500词摘要，足够 S 级分析 |
| bioRxiv 全文 | — | ❌ 多层封锁 | — | JATS XML 403 + PDF 403 + HTML JS渲染 |
| medRxiv API | ✅ 代理 | — | JSON 摘要 | 同 bioRxiv |
| arXiv HTML | ✅ 代理 HTTP-L | ✅ curl /html/{id} | 直取 | 唯一可靠全文源 |
| arXiv API | ✅ 代理 HTTP-L | — | XML 摘要 | — |
| Cell/Nature/Science | ✅ PubMed | ❌ 付费墙 | PubMed摘要 | — |
| Lancet/JAMA | ✅ PubMed | ❌ 付费墙 | PubMed摘要 | — |
| Semantic Scholar | ⚠️ 429 | — | — | 需 API key |

**铁律**: 每次分析论文前，先尝试全文获取（arXiv HTML→可靠；bioRxiv/medRxiv JATS XML→已被 403 封锁）。API 摘要（300-500词）+ LLM 深度推理已足够支撑完整 S 级 7 层解剖——LLM 是基于已有知识体系推理、连接、质疑，不是 NLP 提取。

## Clash 代理

- **端口**: `127.0.0.1:7890` (mixed port)
- **节点**: Ghelper → 美国西雅图
- **规则**: PubMed 直连；arXiv/bioRxiv/medRxiv 走代理
