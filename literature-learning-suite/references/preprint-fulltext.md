# 预印本全文获取技术手册

> 此文档属于 `paper-research` 技能。
> 记录 bioRxiv/medRxiv/arXiv 三源全文获取的实测方案。

## 连通性矩阵 (2026-06-08 更新)

| 源 | Cloudflare | 全文方式 | Cookie | 速度 | 实际可用性 |
|------|:--:|------|:--:|------|:--:|
| **arXiv** | ❌ 无 | `curl /html/{id}` | 不需要 | ~2s/篇 | ✅ 可靠 |
| **bioRxiv API** | ❌ 无 | `api.biorxiv.org/details/...` (JSON摘要) | 不需要 | ~1s/篇 | ✅ 可靠（300-500词摘要） |
| **bioRxiv HTML** | ✅ 有 | `curl` 页面 → **正文JS动态渲染** | `cf_clearance` | ~5s/篇 | ⚠️ CF可过但正文不可提取 |
| **bioRxiv JATS XML** | ✅ 有 | `.source.xml` | `cf_clearance` | — | ❌ 直接返回403（2026-06-08实测） |
| **bioRxiv PDF** | ✅ 有 | `.full.pdf` | `cf_clearance` | — | ❌ 直接返回403（2026-06-08实测） |
| **medRxiv** | ✅ 有 | 同 bioRxiv | `cf_clearance` | ~10s/篇 | ⚠️ 同上限制 |

**关键结论（2026-06-08）**: bioRxiv 全文获取存在**三层壁垒**：
1. **Cloudflare JS Challenge** — `cf_clearance` cookie + 浏览器头可绕过
2. **正文 JS 动态渲染** — 页面是 React/Angular SPA，静态 curl 拿到的是空壳 HTML（~29KB 但无实质文本）
3. **JATS XML / PDF 直接端点 403** — 禁用直接文件访问

**唯一可靠路径**: Playwright MCP（真实 Chrome + JS 执行）或手动浏览器保存。

## 备用方案：Playwright MCP（浏览器自动化，无需 cookie）

**这是目前唯一可靠的 bioRxiv/medRxiv 全文获取方式。**

已安装 `@playwright/mcp` 至 Hermes MCP 配置（`config.yaml` → `mcp_servers.playwright`）。
重启 Hermes 后可用真实 Chrome 浏览器自动通过 Cloudflare **并执行 JS 渲染页面正文**——curl 拿不到 JS 动态加载的内容。

```yaml
# config.yaml
mcp_servers:
  playwright:
    command: npx
    args:
    - "@playwright/mcp@latest"
    - "--browser"
    - "chrome"
    enabled: true
```

启用后流程：`browser_navigate(url)` → Chrome 自动通过 Cloudflare JS 验证 → 页面 JS 渲染正文 → `browser_snapshot(full=true)` 提取全文。

**为何 curl + cookie 不够（2026-06-08 实测）**：
- 即使 `cf_clearance` 有效，bioRxiv 页面是 React SPA，正文通过 JS 动态加载
- curl 拿到的 HTML 只有 ~5KB 骨架（或间歇性 ~29KB 含导航但不含正文）
- `document.querySelector('.article.fulltext')` 在浏览器中才可访问
- JATS XML (`.source.xml`) 和 PDF (`.full.pdf`) 端点直接返回 403

## Chrome Cookie DB 无法程序化读取

Windows 下 Chrome Cookie 存储于 `%LOCALAPPDATA%\Google\Chrome\User Data\Default\Network\Cookies`，
使用 **DPAPI 加密**，无法跨进程解密。不能用 `python -c "import sqlite3; ..."` 读取。

## Cookie 获取（一次配置，数月有效）

bioRxiv 和 medRxiv 共享 Cold Spring Harbor Lab 的 Cloudflare zone。一个 `cf_clearance` 通吃两站。

1. Chrome 打开 `https://www.biorxiv.org/` → 等页面正常加载（自动通过验证）
2. F12 → Application → Cookies → `www.biorxiv.org`
3. 复制 `cf_clearance` 的 Value
4. 如仍使用旧版 cookie 流程，粘贴到 `$LITERATURE_KG_ROOT/.biorxiv_cookie`（纯文本一行）

## curl 必须模拟真实浏览器

仅 `User-Agent` 不够——Cloudflare 还检查 `Sec-Ch-Ua` 和 `Accept`：

```bash
curl -sL "<URL>" \
  -H "Cookie: cf_clearance=XXX" \
  -H "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 ... Chrome/135.0.0.0 Safari/537.36" \
  -H "Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8" \
  -H "Sec-Ch-Ua: \"Google Chrome\";v=\"135\"" \
  -H "Sec-Ch-Ua-Platform: \"Windows\""
```

## 已知陷阱

- **bioRxiv JATS URL 双斜杠**: API 返回的 URL 中有 `//2026`，需正则修复
- **medRxiv 部分 DOI 前缀 API 不支持**: `10.64898` 前缀在 medRxiv API 可能报日期错误，回退直接抓 HTML
- **arXiv 无 Cloudflare**: 不需要 cookie，直接 `/html/{id}`，需代理 7890
- **bioRxiv JATS XML 403**: `.source.xml` 端点已被禁用（2026-06-08），不能用 curl 直接获取
- **bioRxiv PDF 403**: `.full.pdf` 端点同样被禁用
- **bioRxiv HTML 是 SPA**: 正文通过 JS 动态渲染，curl 拿不到实质内容，必须用浏览器（Playwright）
- **cookie 间歇失效**: `cf_clearance` 对同一 IP 的后续请求可能重新触发 CF，同一 session 内 curl + cookie 最多用 1-2 次
- **DOI 前缀兼容性**: `10.64898`（新版）和 `10.1101`（旧版）均可用于 API 检索，但直接访问 URL 以 `10.1101` 为准

## 实际可行的日常流程

对于 cron 定时任务（每日常规深度学习）：
1. **API 摘要** → bioRxiv API JSON（无需 cookie，可靠）→ 300-500 词摘要 → LLM S 级推理
2. **偶尔精读** → 可见 Chrome/CDP → `scripts/extract_biorxiv_cdp.mjs` → 保存到 `fulltext_cache/`
3. **Playwright 自动化**（若已启用）→ `browser_navigate` + `browser_snapshot(full=true)`

摘要 + LLM 推理是规模化可行的——17 篇/天的 cron 任务中，API 摘要已足够支撑完整的 7 层解剖（LLM 不是 NLP 提取，是基于已有知识体系进行推理、连接、质疑）。

## Codex 架构 & Hermes CDP 方案

**Codex CLI 如何成功访问 bioRxiv**：通过 Chrome Native Messaging（`com.openai.codexextension`）
直接操控用户的真实 Chrome 浏览器——Cloudflare 看到的是有 Cookie/登录态/浏览历史的真实用户。

**Hermes 等效方案**：`browser.cdp_url = 'http://127.0.0.1:9222'` + Chrome 远程调试。
详见本 skill 的 `SKILL.md` 中“可见 Chrome / CDP 全文提取”流程。

**快速启动**：`scripts/biorxiv_chrome_cdp_launcher.bat`（需桌面环境，关闭所有 Chrome 窗口后运行）。

## 本 skill 的可用脚本

skill 根目录下和 bioRxiv 相关的脚本：

```bash
scripts/download_biorxiv_api.py          # bioRxiv/medRxiv 官方 API 元数据+摘要
scripts/extract_biorxiv_cdp.mjs          # 从已通过验证的可见 Chrome 页面提取全文
scripts/biorxiv_chrome_cdp_launcher.bat  # 启动独立 CDP Chrome
```

API 下载器用法：

```bash
# 按日期批量下载 API 元数据/摘要（自动翻页）
python scripts/download_biorxiv_api.py --from-date 2026-06-01 --to-date 2026-06-01

# 单篇 DOI 下载 API 记录
python scripts/download_biorxiv_api.py --doi 10.64898/2026.05.31.727600

# medRxiv
python scripts/download_biorxiv_api.py --server medrxiv --from-date 2026-06-01
```

输出 JSONL/CSV/Markdown 到 `$LITERATURE_KG_ROOT/biorxiv_api/`，字段包括 title/authors/doi/category/date/abstract/jatsxml 等。该脚本只抓 API 记录，不绕过 Cloudflare 抓 PDF/JATS 全文。

bioRxiv/medRxiv 自动直连下载全文目前仍不可靠。实际可行路径是：
1. 用 bioRxiv API 获取元数据和摘要；
2. 对需要精读的文章，用 `biorxiv_chrome_cdp_launcher.bat` 启动真实 Chrome；
3. 用 `extract_biorxiv_cdp.mjs` 提取并保存到 `fulltext_cache/`。
