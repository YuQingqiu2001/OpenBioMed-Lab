# MCP 与工具路由

## Hermes 中已发现的 MCP

### PubMed

配置名称：`pubmed`

启动要求：

- 命令：`mcp-simple-pubmed`
- 必填环境变量：`PUBMED_EMAIL`
- 可选：`PUBMED_API_KEY`、`PUBMED_TOOL`

能力：

- `search_pubmed(query, max_results=10)`：检索，最多 50 条。
- `get_paper_fulltext(pmid)`：通过 PMC 获取开放全文。
- 资源：`pubmed://{pmid}/abstract`、`pubmed://{pmid}/full_text`。
- Prompt：系统综述检索、PICO 检索、作者检索。

### arXiv

配置名称：`arxiv`

启动参数：

- 命令：`arxiv-mcp-server`
- `--storage-path <path>`：全文 Markdown、索引和监控状态存储位置。

能力：

- `search_papers`
- `download_paper`
- `list_papers`
- `read_paper`
- `get_abstract`
- `semantic_search`
- `reindex`
- `citation_graph`
- `watch_topic`
- `check_alerts`

推荐链路：`search_papers -> download_paper -> read_paper`。服务端自动遵守 arXiv
至少 3 秒的请求间隔；遇到 429/503 时等待，不要循环重试。

### Fetch

配置名称：`fetch`

能力：抓取普通网页和公开文档。可配置 `--proxy-url`。遵守 robots.txt 和站点访问规则。

### Playwright

配置名称：`playwright`

用于动态网页、JavaScript 渲染正文和用户可见浏览器流程。对 bioRxiv/medRxiv
动态正文优先使用它。不得自动破解 CAPTCHA/Turnstile；需要时由用户手动完成。

## 配置模板

见 `assets/templates/mcp-servers.yaml`。替换：

- `${PUBMED_EMAIL}`
- `${ARXIV_STORAGE_PATH}`

不要把真实邮箱、API key、cookie 或代理认证写进可分发 skill。

## 回退原则

1. MCP 可用时优先 MCP。
2. MCP 不可用时使用本 skill 的 CLI 脚本。
3. 网络源失败时记录错误、日期和已尝试路径。
4. 不把 0 条结果自动解释为“没有相关论文”；先排除网络、查询语法和限流问题。
