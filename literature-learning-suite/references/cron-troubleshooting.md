# Cron 定时任务排障指南

## 症状：任务从未运行 (last_run_at: null)

### 排查步骤

#### 1. 检查任务是否存在且启用
在 Hermes 会话中: `cronjob(action='list')`
确认 `enabled: true`, `state: scheduled`。

#### 2. 检查时区（最常见根因 2026-06-07）

**症状**: `next_run_at` 显示 `+01:00` 而非 `+08:00`(HKT)

**根因**: `config.yaml` 中 `timezone: ''`(空) → 调度器默认 UTC+1

**修复**: `config.yaml` → `timezone: 'Asia/Hong_Kong'`

**重启前变通** (实测有效):
| 调度器时区 | 目标 HKT | cron 表达式 |
|-----------|----------|-------------|
| UTC+1 | 09:00 HKT | `0 2 * * *` |
| UTC | 09:00 HKT | `0 1 * * *` |

重启后恢复 `0 9 * * *`。

#### 3. 检查创建时间
若 `created_at` 晚于当天的触发时间,该天被跳过,次日首次运行。

#### 4. 检查调度器进程
Windows: `tasklist | grep -i hermes`
Lock 文件: `~/AppData/Local/hermes/cron/.tick.lock`

#### 5. 强制触发
`cronjob(action='run', job_id='81d87c68a429')`

### 日志
```bash
~/AppData/Local/hermes/logs/agent.log   # 搜索 cron_<job_id>
~/AppData/Local/hermes/logs/errors.log
```

### 已知限制
timezone 变更需重启 Hermes 才被调度器拾取。重启前用 cron 表达式补偿。
