# S-tier 批量升级工作流

> 2026-06-08 实战验证。将 NLP 提取模式的论文升级到 v4.0 真正 S 级 7 层解剖。

## 背景

知识图谱中多数论文使用 NLP 提取模式生成（2026-06-07 之前）：
- **T3**: 1 条 C-E-S 链，claim = evidence = 摘要复制
- **T4**: trigger = 论文标题，cascade = 摘要字符串
- **T5**: observation_1 = 摘要第一句
- **T6**: new_concepts = 论文标题
- **T7**: 1 条关联，全为 `same_journal`（非生物边，已禁止）

升级后每条论文需要：
- T3: 5 条独立 C-E-S 链（各有独立 claim/evidence/synthesis/strength/uncertain）
- T4: 完整机制级联（trigger → cascade → key_modifications → feedback）
- T5: 3 组 observation → interpretation
- T6: new_concepts / overturned_views / methodological_breakthroughs
- T7: 5 条有实质描述的关系（非 same_journal）
- analysis_method: `"LLM_deep_reasoning_S_tier_v4.0"`

## 审计

```bash
# 完整自检
python D:/knowledge_graph/scripts/selfcheck_knowledge_graph.py

# 快速查看 weak 论文分布
python -c "
import json
papers = [json.loads(l) for l in open('./literature-workspace/papers.db')]
s = [p for p in papers if p.get('analysis_tier')=='S']
t3 = sum(1 for p in s if len(p.get('tier3_ces_chains',[]))>=5)
t7 = sum(1 for p in s if len(p.get('tier7_cross_refs',[]))>=5)
print(f'Total S: {len(s)}; T3>=5: {t3}; T7>=5: {t7}; Strict pass: {sum(1 for p in s if len(p.get(\"tier3_ces_chains\",[]))>=5 and len(p.get(\"tier7_cross_refs\",[]))>=5)}')
"
```

## 批量升级模式

### 单篇模板

```python
U['PMID:XXXXXXXX'] = {
    'tier3_ces_chains': [
        {'chain_id':1,'claim':'...','evidence':'... specific data, numbers, methods ...',
         'synthesis':'... integrative interpretation ...','strength':3,
         'uncertain':'... limitations ...'},
        # ... 5 chains total
    ],
    'tier4_mechanism_cascade': {
        'trigger': '...',
        'cascade': ['step1','step2',...],
        'key_modifications': [{'site':'...','mod':'...','effect':'...'}],
        'downstream_effects': '...',
        'feedback': [{'type':'positive','node':'...'}]
    },
    'tier5_hidden_axis': [
        {'observation':'...','interpretation':'...'},
        ... 3 pairs
    ],
    'tier6_concept_innovation': {
        'new_concepts': ['...'],
        'overturned_views': ['...'],
        'methodological_breakthroughs': ['...']
    },
    'tier7_cross_refs': [
        {'ref_id':'PMID:YYYYYYYY','relation':'shared_...','description':'...'},
        ... 5 refs
    ],
    'analysis_method':'LLM_deep_reasoning_S_tier_v4.0',
    'analysis_updated':'2026-06-08',
}
```

### 执行步骤

```bash
# 1. write_file 创建升级脚本
# 2. terminal 执行
python D:/knowledge_graph/scripts/upgrade_batchN.py

# 3. 重建边和网络
rm -rf D:/knowledge_graph/scripts/__pycache__
python -B D:/knowledge_graph/scripts/gen_edges.py
python -B D:/knowledge_graph/scripts/build_network.py

# 4. 验证
python D:/knowledge_graph/scripts/selfcheck_knowledge_graph.py
```

## ⚠️ 关键陷阱

### Windows MSYS heredoc

**症状**: `terminal` 中使用 `python - <<'PY'` 长脚本时，出现 `unexpected EOF`。

**原因**: MSYS bash 对含嵌套引号的 Python heredoc 处理错误。

**解决**: 永远用 `write_file` → `terminal("python script.py")` 两段式，不在终端直接写长 Python。

### gen_edges.py tier5 列表格式

**症状**: `AttributeError: 'list' object has no attribute 'values'`

**原因**: gen_edges.py 策略 4（line ~424）假设 `tier5_hidden_axis` 是 dict 格式 `{observation_1, interpretation_1}`，但 v4.0 升级使用 list 格式 `[{observation, interpretation}]`。

**修复**: gen_edges.py v3.1+ 已修复（添加 `_ha_text()` 兼容函数），检查 line ~424-433。

### T7 跨文献关联的边生成

升级后的 T7 使用 `ref_id`/`relation`/`description` 格式（非旧 `pmid`/`relation`/`desc`），gen_edges.py 策略 1 自动提取这些关系生成新边。

## 实战结果

| 批次 | 论文数 | 升级前 s_weak | 升级后 s_weak | 新增边 |
|------|--------|--------------|-------------|--------|
| Batch 1 | 10 | 449 | 439 | 107 |
| Batch 2 | 10 | 439 | 429 | 154 |
| **累计** | **20** | **449** | **429** | **261** |

剩余 429 篇待升级。按每批 10 篇计算，约需 43 批次完成全部 449 篇的 v4.0 升级。
