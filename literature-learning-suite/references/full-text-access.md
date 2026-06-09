# 文献全文获取策略

所有命令默认从 skill 根目录运行，并先设置：

```powershell
$env:LITERATURE_KG_ROOT = "D:\path\to\literature-workspace"
```

## 1. arXiv

先检索元数据：

```bash
python scripts/search_arxiv.py "single cell foundation model" --max-results 20
```

下载 PDF 后优先使用 PyMuPDF 提取：

```bash
python scripts/extract_pymupdf.py paper.pdf --output "$LITERATURE_KG_ROOT/fulltext_cache/ARXIV_2402.03300.txt"
```

复杂版面、公式或多栏顺序异常时使用 Marker：

```bash
python scripts/extract_marker.py paper.pdf --output-dir "$LITERATURE_KG_ROOT/fulltext/marker"
```

## 2. PubMed / PMC

通过 NCBI E-utilities 将 PMID 映射为 PMCID，再获取开放全文。只有进入 PMC 的开放获取论文才保证有全文；其余记录应保留摘要与来源链接，不要伪造全文。

```text
PMID -> elink(dbfrom=pubmed, db=pmc) -> PMCID -> PMC XML/HTML
```

将清洗后的正文写入：

```text
$LITERATURE_KG_ROOT/fulltext_cache/PMID_42248141.txt
```

## 3. bioRxiv / medRxiv 元数据与摘要

官方 API 适合发现、筛选和摘要分析：

```bash
python scripts/download_biorxiv_api.py --from-date 2026-06-01 --to-date 2026-06-01
python scripts/download_biorxiv_api.py --doi 10.64898/2026.05.31.727600
python scripts/download_biorxiv_api.py --server medrxiv --from-date 2026-06-01
```

输出位于 `$LITERATURE_KG_ROOT/biorxiv_api/`。API 记录不是全文；分析时必须明确区分摘要证据和全文证据。

## 4. bioRxiv 可见 Chrome / CDP

当页面需要 JavaScript 渲染或 Cloudflare 验证时，使用已经在 Hermes 本机流程中验证过的可见 Chrome 路径。

1. 关闭占用目标 profile 的 Chrome。
2. 运行 `scripts/biorxiv_chrome_cdp_launcher.bat`。
3. 在可见窗口中完成必要验证，但不要登录 Google。
4. 提取单篇或批量全文：

```bash
node scripts/extract_biorxiv_cdp.mjs --doi 10.64898/2026.05.31.727600 --port 9223
node scripts/extract_biorxiv_cdp.mjs --batch doi_list.txt --port 9223
```

独立 profile 默认位于 workspace 内。它可能保存 Cookie；不要提交 profile、Cookie 或浏览器状态到 Git。

## 5. 手动回退

如果自动提取仍失败：

1. 在可见浏览器中打开并等待正文完整渲染。
2. 保存页面或复制正文。
3. 清洗导航、参考文献噪声和重复段落。
4. 写入 `fulltext_cache/`，并保留来源 URL、抓取日期和提取方法。

## 缓存命名

```text
$LITERATURE_KG_ROOT/fulltext_cache/{normalized_id}.txt
```

规范化规则：

```python
safe_id = re.sub(r"[^a-zA-Z0-9_.-]", "_", paper_id)
```

示例：

```text
PMID:42248141               -> PMID_42248141.txt
ARXIV:2402.03300            -> ARXIV_2402.03300.txt
BIORXIV:10.1101/2026.01.01  -> BIORXIV_10.1101_2026.01.01.txt
```

`scripts/gen_edges.py` 会优先读取该缓存来计算语义与实体关联。
