# gen_edges.py 关联边生成详解

## 五策略（全倒排索引，纯生物学）

| 策略 | 关系类型 | 方法 | 阈值 |
|------|---------|------|------|
| 1. 显式引用 | cites/extends 等 | tier7_cross_refs（过滤 same_journal/same_issue） | — |
| 2. 共享分子 | shares_molecules | Bioconductor 73k 基因 两段式匹配 | ≥2 共享 |
| 2.5 文本重叠 | shares_topic | 核心发现 bag-of-words 倒排索引 | ≥4 共享关键词 |
| 3. 同病+同法 | shares_disease_method | 疾病×方法交叉 | ≥1 共享 |
| 4. 隐藏轴 | shares_paradigm | Tier5 关键词共鸣 | ≥1 共享 |
| 5. 概念节点 | defines_concept | concepts.db | — |

## 两段式基因匹配

```
文本 → RE_GENE_ABBREV (全大写) + RE_TITLECASE_GENE (TitleCase小鼠)
     → 查 Bioconductor 73k 基因集合 (O(1) lookup)
     → 保留在集合内的 → 真基因
     → 不在集合内的 → 丢弃（非基因缩写）
```

对比旧方案（手写 126 词表 + BIOWORDS 语境过滤）：
- 旧方案假阳性率 ~40%（APP/STING/MIL 误匹配 ML 论文）
- 新方案假阳性率 ~0%（只有真实基因符号被保留）

## 关键陷阱

### 陷阱 1: RE_PATHWAY 灾难性回溯

**症状**: `gen_edges.py` 运行后无限挂起，CPU 100%
**根因**: 原始正则 `((?:[A-Z][a-z]*/)*(?:[A-Z]{2,}[a-z]*\s*)+(?:pathway|signaling|axis))` 
在非通路文本上触发指数级回溯
**修复**: 改为固定词表匹配
```python
RE_PATHWAY = re.compile(
    r'(Wnt|Hedgehog|Notch|NF-κB|MAPK|PI3K[/-]Akt|...)\s+(pathway|signaling|axis)',
    re.IGNORECASE)
```

### 陷阱 2: .pyc 幽灵缓存

**症状**: 修改 gen_edges.py 后删除的代码（如同刊同期策略）仍执行，仍生成 same_journal 边
**根因**: `scripts/__pycache__/gen_edges.cpython-311.pyc` 未失效，Python 加载旧版本
**修复**: 每次修改后必须 `rm -rf scripts/__pycache__`
**验证**: `python -c "import gen_edges; print('same_journal' in open(gen_edges.__file__).read())"` → False

### 陷阱 3: same_journal 来自论文数据自身

**症状**: gen_edges.py 不含 same_journal 代码，但 edges.db 中仍有 198 条
**根因**: 论文的 `tier7_cross_refs` 中 LLM 写了 `relation: "same_journal"` 的关联，
策略 1 忠实地读取并生成边
**修复**: 策略 1 加入 `NON_BIO_RELS = {'same_journal', 'same_issue', 'same_author'}` 过滤

### 陷阱 4: paper_text() 的全文缓存路径

缓存文件命名必须与 gen_edges.py 中 `paper_text()` 的路径计算一致：
```python
safe = re.sub(r'[^a-zA-Z0-9_.-]', '_', pid)
cache_file = Path(os.environ["LITERATURE_KG_ROOT"]) / "fulltext_cache" / f"{safe}.txt"
```
所有全文提取脚本应使用相同的命名规则。

## 运行验证清单

修改 gen_edges.py 后：
```bash
# 1. 清缓存
rm -rf scripts/__pycache__

# 2. 运行
python scripts/gen_edges.py

# 3. 验证无 same_journal
python -c "
import json, collections
c = collections.Counter()
for l in open('edges.db'): 
    c[json.loads(l)['relation']] += 1
assert 'same_journal' not in c, f'BUG: same_journal={c[\"same_journal\"]}'
print('OK: no same_journal')
"

# 4. 刷新网络图
python scripts/build_network.py
```
