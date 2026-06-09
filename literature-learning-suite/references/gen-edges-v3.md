# gen_edges.py v3.1 — 生物学关联边生成架构

> 关联 gen_edges.py v3.1, 2026-06-07 | 取代手写 126 词表，接入 Bioconductor 73k 基因

## 架构总览

```
输入: 492 篇论文 + 11 概念 + 73,017 基因符号 + 297 KEGG 通路

策略 1  显式引用  ──→ tier7 手工关联（排除 same_journal/same_issue/same_author）
策略 2  共享分子  ──→ 两篇论文共享 ≥2 基因/通路（Bioconductor 验证）
策略 2.5 文本重叠 ──→ 核心发现关键词重合 ≥4（主力，~400 条）
策略 3  同病+同法 ──→ 疾病 × 技术 交叉（倒排索引）
策略 4  隐藏轴   ──→ Tier5 深层范式共鸣（paradigm/bias/survivor 等关键词）
策略 5  概念节点  ──→ 论文定义 → 概念

输出: ~480 条纯生物学边
```

## 倒排索引 > O(N²) 全对组合

所有策略统一使用倒排索引（entity → paper_ids），避免 492² = 242k 次比较：

```python
# 策略 2 示例: 基因 → 论文
mol_index = defaultdict(set)           # "TNF" → {"PMID:1","PMID:2","PMID:3"}
for p in deep_papers:
    for gene in extract_entities(p):   # Bioconductor 两段式匹配
        mol_index[gene].add(p.id)

# 只需为出现在同一基因下的论文建对（而非所有组合）
for gene, pids in mol_index.items():
    for p1, p2 in combinations(pids, 2):  # 每基因最多 15 篇
        if shared_count >= 2: add_edge(...)
```

## Bioconductor 基因集导入

### 数据源
- `org.Hs.eg.db` → 33,229 人类基因符号（全大写：TNF, EGFR, BCL2）
- `org.Mm.eg.db` → 56,913 小鼠基因符号（TitleCase: Tnf, Il6, Cd4, Trp53）
- 合并去重 → **73,017 个唯一基因符号**

### 导出命令（WSL）
```bash
wsl -d Ubuntu-20.04 -e bash -c "
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate sc_spatial_env
Rscript -e '
library(org.Hs.eg.db); library(org.Mm.eg.db); library(jsonlite)
hs <- grep(\"^[A-Z][A-Z0-9]{1,7}$\", keys(org.Hs.eg.db, \"SYMBOL\"), value=TRUE)
mm <- grep(\"^[A-Z][a-z0-9]{1,7}$\", keys(org.Mm.eg.db, \"SYMBOL\"), value=TRUE)
write_json(unique(c(hs, mm)), "bioc_genes.json")
'
"
```

### 两段式匹配（O(候选数)，非 O(词表×文本)）
```python
# 第一段: 正则提取候选
RE_GENE_ABBREV = re.compile(r'\b([A-Z]{2,6}\d{0,3}[A-Z]?)\b')     # TNF, EGFR, BCL2
RE_TITLECASE_GENE = re.compile(r'\b([A-Z][a-z]{2,7}\d{0,3}[a-z]?)\b')  # Tnf, Il6, Smad2

# 第二段: 查 73k 基因集合 O(1)
genes = load_json("bioc_genes.json")
for candidate in extract_candidates(text):
    if candidate.upper() in genes:   # 是真基因 → 保留
        found.append(candidate)
```

## ⚠️ 关键陷阱

### 陷阱 1: `.pyc` 缓存幽灵
**症状**: gen_edges.py 修改后仍然生成 `same_journal` 边，尽管源码已删除相关代码。

**根因**: Python 的 `__pycache__/gen_edges.cpython-311.pyc` 保留了旧版本字节码。新进程 import 时优先使用 .pyc，忽略源码修改。

**修复**:
```bash
rm -rf D:/knowledge_graph/scripts/__pycache__
python -B scripts/gen_edges.py    # -B 禁止写入 .pyc
```

### 陷阱 2: 灾难性正则回溯
**症状**: `RE_PATHWAY` 正则在非通路文本上 O(e^n) 回溯，单篇论文耗时 >30s。

**问题正则**:
```python
# ❌ 嵌套量词导致回溯爆炸
re.compile(r'\b((?:[A-Z][a-z]*/)*(?:[A-Z]{2,}[a-z]*\s*)+(?:pathway|signaling))\b')
```

**修复**: 改为固定词表匹配
```python
# ✅ 安全: 具体词表 OR
re.compile(r'(Wnt|Hedgehog|Notch|MAPK|PI3K[/-]Akt|mTOR|AMPK|...)\s+(pathway|signaling)')
```

### 陷阱 3: tier7_cross_refs 含非生物边
**症状**: 从空 edges.db 重建，仍出现 198 条 `same_journal` 边。

**根因**: 策略 1 读取 `papers.db` 中每篇论文的 `tier7_cross_refs` 字段——LLM 在 S 级分析时写了 `relation: "same_journal"` 的关联。gen_edges.py 策略 1 忠实地把它们作为显式引用边生成。

**修复**: 策略 1 加黑名单过滤
```python
NON_BIO_RELS = {'same_journal', 'same_issue', 'same_author'}
if ref.get('relation', 'cites') in NON_BIO_RELS:
    continue
```

### 陷阱 4: 高频基因配对爆炸
**症状**: TNF/NF-κB 出现在 100+ 篇论文中 → 倒排索引生成 100²/2 = 4,950 对。

**修复**: 每基因最多取 15 篇
```python
MAX_PAPERS_PER_MOL = 15
pid_list = list(pids)[:MAX_PAPERS_PER_MOL]
```

## KEGG 通路过滤

KEGGREST 导出 371 条通路名称中混入疾病通路（"Colorectal cancer", "Pancreatic cancer" 等）→ 文本匹配造成误识别。

过滤规则：
```python
disease_terms = ['cancer','carcinoma','tumor','melanoma','tuberculosis',
                 'leukemia','lymphoma','diabetes','alzheimer','parkinson',
                 'infection','disease','virus','hepatitis','influenza','...']
clean = [p for p in pathways if not any(dt in p.lower() for dt in disease_terms)]
# 371 → 297
```

## 关联边类型一览

| 关系类型 | 来源策略 | 说明 |
|----------|----------|------|
| `shares_topic` | 2.5 | 核心发现关键词重叠 ≥4 |
| `shares_molecules` | 2 | Bioconductor 基因共享 ≥2 |
| `shares_focus` | 1 | 论文数据中标签共享 |
| `shares_paradigm` | 4 | Tier5 隐藏轴共鸣 |
| `shares_disease_method` | 3 | 同疾病 + 同技术 |
| `extends` | 1 | 机制延伸（如 FIDELIO→FIGARO 汇总分析） |
| `accompanied_by` | 1 | 述评配主文 |
| `defines_concept` | 5 | 概念定义 |
| `same_class` | 1 | 同类药理 |

**已删除**: `same_journal`(198), `same_issue`(62) — 非生物学关联。
