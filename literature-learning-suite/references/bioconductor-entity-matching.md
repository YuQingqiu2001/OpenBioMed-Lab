# Bioconductor 实体匹配 — gen_edges.py v3.1 核心技术

## 问题

v3.0 手写 126 个分子实体（KNOWN_MOLECULES）覆盖率不足：
- 只能识别 ~30% 的论文间分子关联
- 动态大写缩写匹配（RE_GENE_ABBREV）产生大量假阳性：APP(Application) / STING(动词) / MIL(Multiple Instance Learning) 被误当基因
- BIOWORDS 语境过滤器可减少但无法根除假阳性

## 方案：Bioconductor 词表 + 两段式匹配

### 架构

```
org.Hs.eg.db (R/Bioconductor)  ──R export──►  bioc_genes.json (33k 基因)
KEGGREST (R)                    ──R export──►  kegg_pathways.json (297 通路)
                                                    │
                                                    ▼
论文文本  ──RE_GENE_ABBREV──►  候选缩写  ──查集合──►  真实基因 ✓
                                                    │
                                         不在集合中 ──►  非基因 ✗ (丢弃)
```

### 两段式匹配原理

```
# 第一段: 提取候选 (快速正则)
candidates = RE_GENE_ABBREV.findall(text)  # e.g., ['TNF', 'APP', 'EGFR', 'MIL']

# 第二段: 查权威集合 (O(1) per candidate)
genes = load_bioc_genes()  # set of 33,229 gene symbols
valid_genes = {c for c in candidates if c in genes}
# → {'TNF', 'EGFR'} — APP 和 MIL 不在基因集合中，自动滤除
```

### 复杂度

| 方案 | 时间复杂度 | 本次实测 |
|------|-----------|---------|
| 手写遍历 (v3.0) | O(词表 × 文本) | 0.1s (仅126词) |
| **两段式 (v3.1)** | O(候选数) + O(set加载) | **0.7s** (33k词) |

## 灾难性回溯修复

### 问题

v3.0 的 `RE_PATHWAY` 正则使用了嵌套量词：

```python
# ❌ 灾难性回溯 — 0.7s → 永不完成
RE_PATHWAY = re.compile(
    r'\b((?:[A-Z][a-z]*/)*(?:[A-Z]{2,}[a-z]*\s*)+'
    r'(?:pathway|signaling|pathway|axis))\b', re.IGNORECASE
)
```

在 "Wnt/planar cell polarity pathway activation triggers..." 这样的文本上，
正则引擎尝试数百万种解析路径后超时。

### 修复

```python
# ✅ 固定词表 — O(1) per match
RE_PATHWAY = re.compile(
    r'(Wnt|Hedgehog|Notch|TGF-β|NF-κB|MAPK|PI3K[/-]Akt|JAK[/-]STAT|mTOR|AMPK|'
    r'Hippo|RAS|RAF|RTK|GPCR|JNK|p38|ERK|cGAS-STING|inflammasome|complement|'
    r'ferroptosis|autophagy|apoptosis|pyroptosis|necroptosis|'
    r'glycolysis|OXPHOS|PPP|TCA)\s+(pathway|signaling|axis|signalling|cascade)',
    re.IGNORECASE
)
```

### 排查方法

```python
import time
# 分段计时，定位慢的正则
for name, pat in [("RE_GENE_ABBREV", RE_GENE_ABBREV), 
                   ("RE_PATHWAY", RE_PATHWAY)]:
    t0 = time.time()
    pat.findall(test_text)
    print(f"{name}: {time.time()-t0:.3f}s")
```

RE_GENE_ABBREV < 0.001s → RE_PATHWAY 超时 → 问题定位。

## KEGG 通路过滤

KEGG 包含疾病通路（"Colorectal cancer", "Pancreatic cancer"），
在论文文本中匹配这些疾病名会产生虚假分子关联。

### 过滤规则

```python
disease_terms = ['cancer','carcinoma','tumor','melanoma','tuberculosis',
                 'leukemia','lymphoma','infection','disease','virus', ...]
clean_paths = [p for p in paths 
               if not any(dt in p.lower() for dt in disease_terms)]
# 371 → 297 条
```

## 倒排索引优化

全对组合 O(N²) 在 484 篇论文下超时（484² × 分子匹配 = 数百万次比较）。

### 修复: 分子→论文 倒排索引

```python
# ❌ O(N²)
for p1, p2 in combinations(papers, 2):
    shared = p1.molecules & p2.molecules
    if len(shared) >= 2: add_edge(p1, p2)

# ✅ 倒排索引 O(M) where M = number of unique molecules
mol_index = defaultdict(set)  # mol → {paper_ids}
for p in papers:
    for mol in p.molecules:
        mol_index[mol].add(p.id)

# 只在共享同一分子的论文间检查
pair_shared = defaultdict(set)
for mol, pids in mol_index.items():
    for p1, p2 in combinations(list(pids)[:15], 2):  # 每分子限15篇
        pair_shared[(p1,p2)].add(mol)
```

加上每分子限 15 篇（MAX_PAPERS_PER_MOL=15），避免 TNF/NF-κB 等高频分子产生百万级组合。

## 小鼠基因 TitleCase 匹配 (v3.1 新增)\n\n### 问题\n\n人类基因命名规范为全大写（`TNF`, `EGFR`, `BCL2`），但**小鼠基因**使用 TitleCase（`Tnf`, `Egfr`, `Bcl2`）。\n原 `RE_GENE_ABBREV = r'\\b([A-Z]{2,6}\\d{0,3}[A-Z]?)\\b'` 只能匹配全大写，完全遗漏了小鼠基因。\n\n### 修复\n\n新增 `RE_TITLECASE_GENE` 正则，添加 `org.Mm.eg.db` 基因到词表：\n\n```python\n# 全大写 — 人 + 部分小鼠 (C2, C3, F5, etc.)\nRE_GENE_ABBREV = re.compile(r'\\b([A-Z]{2,6}\\d{0,3}[A-Z]?)\\b')\n\n# TitleCase — 小鼠主流 (Tnf, Il6, Cd4, Trp53, Smad2, Rac1)\nRE_TITLECASE_GENE = re.compile(r'\\b([A-Z][a-z]{2,7}\\d{0,3}[a-z]?)\\b')\n```\n\n两条正则提取的候选统一通过 `upper()` 后查 73k 基因集合。TitleCase 模式需 ≥3 个小写字母\n（避免匹配 `The`/`This`/`And` 等常见英语单词），即使误匹配也会被基因集合滤除。\n\n### 词表规模（人+鼠合并）\n\n| 来源 | 数据库 | 基因数 | 模式 |\n|------|--------|--------|------|\n| 人类 | `org.Hs.eg.db` | 33,229 | ALLCAPS |\n| 小鼠 | `org.Mm.eg.db` | 56,913 | TitleCase + 少量 ALLCAPS |\n| **合并去重** | | **73,017** | |\n\n### TitleCase 发现的 20 个新基因\n\nSMAD2, RAC1, GDF15, PIEZO2, GAS6, ARG1, FOS, IRF4, RAB7, SKP2, CRK, DES, FER, AXL, 等\n\n产生了高质量关联边（例：SMAD2+TGF-β → FSTL1⇄IGFBP4 共享 TGF-β 信号轴）。\n\n## 陷阱

### .pyc 缓存幽灵执行

修改 `gen_edges.py` 后删除 `scripts/__pycache__/`，否则旧 `.pyc` 会继续执行已删除的代码。
**实测**: 删除了同刊同期生成代码后，`same_journal` 边仍然出现——最终定位到 `.pyc` 未失效。
`python -B` 可临时跳过 .pyc，但常规用法是 `rm -rf scripts/__pycache__`。

### 僵尸进程

`terminal()` 超时时子进程不会被杀。若正则陷入灾难性回溯，进程持续消耗 CPU。
排查: `ps aux | grep python` → `kill <PID>`。
本次清理: 4 个 `gen_edges.py` 僵尸 (PID 11314/11353/11473/11481) 来自 RE_PATHWAY 回溯。

```bash
# 1. WSL 中运行 R 导出
wsl -d Ubuntu-20.04 -e bash -c "
  source /path/to/miniconda3/etc/profile.d/conda.sh
  conda activate sc_spatial_env
  Rscript D:/knowledge_graph/scripts/export_bioc_genes.R
"

# 2. 重新生成关联边
cd D:/knowledge_graph && python scripts/gen_edges.py

# 3. 刷新网络图
cd D:/knowledge_graph && python scripts/build_network.py
```
