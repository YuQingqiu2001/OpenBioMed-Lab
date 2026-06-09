# 代码自检清单 — 修改 KG 脚本后必查

## 流水线端到端验证

每次修改 `gen_edges.py` / `build_network.py` / `kg_core.py` / 全文提取脚本后执行：

```bash
cd ./literature-workspace

# 1. 清除 .pyc 缓存（关键！否则旧代码可能仍执行）
rm -rf scripts/__pycache__

# 2. gen_edges
python scripts/gen_edges.py
# 检查: 无 Exception, 输出 "已有边: N + 新增: M = 总计: T"
# 检查: 新增边类型中无 same_journal / same_issue

# 3. build_network
python scripts/build_network.py
# 检查: "✅ ./literature-workspace/network.html (XXX KB)"


# 4. gen_digest
python scripts/gen_digest.py
# 检查: "日报已生成: ./literature-workspace/daily_digest/YYYY-MM-DD.md"

# 5. kg_core 完整性
python -c "from scripts.kg_core import *; s=get_stats(); print(s)"
# 检查: papers / concepts / edges 数字正确
```

## gen_edges.py 专项检查

| 检查项 | 命令 | 通过标准 |
|--------|------|----------|
| .pyc 已清除 | `ls scripts/__pycache__/` 2>/dev/null | 空或无此目录 |
| same_journal 不存在 | `python -c "import json; rels=[json.loads(l)['relation'] for l in open('./literature-workspace/edges.db') if l.strip()]; print('same_journal' in rels)"` | `False` |
| 策略覆盖 | 看 gen_edges 输出 | 策略1-5 均有输出行 |
| 分子实体覆盖 | `python -c "..." ` | >100 unique molecules across papers |
| 通路匹配可用 | `ls -lh data/kegg_pathways.json` | >100KB (含 25k+ GO 术语) |

## .pyc 缓存陷阱（2026-06-07 实测）

**症状**: 删除了 `gen_edges.py` 中生成 `same_journal` 边的代码后，重新运行仍生成 198 条 `same_journal` 边。

**根因**: `scripts/__pycache__/gen_edges.cpython-311.pyc` 包含旧代码。Python 在源文件修改时间未变化时优先使用 .pyc。

**修复**: `rm -rf scripts/__pycache__` 后重跑。为保险，也可 `python -B scripts/gen_edges.py`（-B 禁止写 .pyc）。

## same_journal 幽灵边的双重根因（2026-06-07 深入排查）

**症状升级**: 清除 .pyc 后仍然生成 same_journal 边。

**第二根因**: 不是 gen_edges.py 的策略代码生成的——是 **papers.db 中 tier7_cross_refs 字段**已包含 LLM 写的 `relation: "same_journal"` 关联。策略1（显式引用）忠实地读取出这些边。

**修复**: 策略1 增加 `NON_BIO_RELS = {'same_journal', 'same_issue', 'same_author'}` 过滤。

## kg_core.py 历史 bug 清单

| 时间 | bug | 修复 |
|------|-----|------|
| 2026-06-07 | `ensure_dirs` 创建 `papers/`/`concepts/`/`relationships/` 空目录 | 改为仅创建 `KG_ROOT` + `DIGEST_DIR` |
| 2026-06-07 | `generate_paper_id` 缺失 `arxiv_id` → `ARXIV:` 前缀 | 增加 `if paper.get("arxiv_id"): return f"ARXIV:{...}"` |
| 2026-06-07 | `lookup_journal_impact_factor` 模糊匹配 "Cell" 误匹配 "Cell Reports" | 增加长度比 >0.7 过滤: `ratio = min(len(norm),len(key)) / max(len(norm),len(key))` |
| 2026-06-07 | `gen_digest.py` 硬编码日期 `TODAY = "2026-06-07"` + 硬编码统计 | 重写为动态读取 `papers.db`/`edges.db` |

## 常见静默失败

| 症状 | 可能原因 |
|------|----------|
| 边数为 0 | edges.db 格式错误（如写入了 `[]`），清理或重建 |
| 新增边全为 shares_topic | 策略2/3/4 的输入为空——检查 deep_papers 过滤条件 |
| 指纹完成 491 而非 492 | 一篇论文 JSON 解析失败，检查 `papers.db` 最后一行 |
| 网络图节点数异常 | edges.db 中 source/target 引用的 ID 不在 papers.db 中 |
| `read_file` 拒绝 papers.db | .db 扩展名被识别为二进制，用 `terminal` + `python -c` 替代 |

## cron 任务验证

```bash
# 手动触发（返回 success=true 仅表示入队成功，不代表已执行完成）
# 等5-10分钟后检查：
wc -l ./literature-workspace/papers.db  # 看有无新增
cronjob(action='list')               # 看 last_status
```

`last_status: "ok"` 且 `papers.db` 行数增加 = 成功。`last_delivery_error` 是投递层问题，不影响数据持久化。
