# S级自检流程

> 强制执行——每次写入 papers.db 后运行，验证非空壳。

## 快速检查

```bash
cd D:/knowledge_graph && python -c "
import json
p=[json.loads(l) for l in open('papers.db','r') if l.strip()]
today=[x for x in p if x.get('analysis_date')=='2026-06-07']

# Tier 3 evidence check
t3_ok=sum(1 for x in today if x.get('tier3_ces_chains') and len(x['tier3_ces_chains'])>0 and x['tier3_ces_chains'][0].get('evidence','') and len(x['tier3_ces_chains'][0]['evidence'])>20)
t3_total=sum(1 for x in today if x.get('tier3_ces_chains'))

# Tier 4 mechanism check
t4_ok=sum(1 for x in today if x.get('tier4_mechanism_cascade',{}).get('cascade') and len(x['tier4_mechanism_cascade']['cascade'])>0)

print(f'今日: {len(today)}篇')
print(f'T3有实质evidence: {t3_ok}/{t3_total}')
print(f'T4有机制级联: {t4_ok}/{len(today)}')

if t3_ok < t3_total:
    print('⚠️ 有空壳T3——检查 tier3_ces_chains[0].evidence 字段')
if t4_ok < len(today):
    print('⚠️ 有空壳T4——检查 tier4_mechanism_cascade 字段')
"
```

## 空壳检测标准

| 字段 | 空壳标志 | 修复 |
|------|----------|------|
| `tier3_ces_chains[0].evidence` | 缺失或 <20 字符 | 用 `core_findings` 填充 evidence |
| `tier3_ces_chains[0].synthesis` | 缺失 | 结合 claim + findings 生成 |
| `tier4_mechanism_cascade.cascade` | 缺失或空数组 | 用 `mechanism_cascade` 字符串填充 |
| `tier5_hidden_axis` | 缺失 | 用 `core_findings` 作为 observation_1 |
| `tier7_cross_refs` | 空数组 | 基于共享标签/期刊自动生成 |

## 严禁事项

- ❌ 只改 `analysis_tier: "S"` 标签而不填充 Tier 2-7 字段
- ❌ 批量套壳后不检查就宣布"100%完成"
- ❌ 用 `sq1: "待深入"` / `sq1: "detail"` 等占位符填充子问题
- ✅ 每次写入后立即运行自检脚本
