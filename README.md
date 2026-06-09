# OpenBioMed Lab

> 开放生物医学人工智能研究实验室
>
> 可复现、可分发、可协作的 AI 驱动生物医学研究工具集合。
> 每一行代码都经过真实研究场景验证。

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-blue.svg)](LICENSE)
[![Hermes Agent](https://img.shields.io/badge/Powered_by-Hermes_Agent-8b5cf6)](https://github.com/NousResearch/hermes-agent)
[![DeepSeek](https://img.shields.io/badge/LLM-DeepSeek_V4_Pro-4f46e5)](https://www.deepseek.com/)

---

## 目录

- [实验室简介](#实验室简介)
- [推荐技术栈](#推荐技术栈)
- [子项目](#子项目)
- [安装与使用](#安装与使用)
- [许可证](#许可证)
- [路线图](#路线图)

---

## 实验室简介

**OpenBioMed Lab** 是一个开放的生物医学 AI 研究工具仓库。这里的每一个子项目
都源自真实的医学 AI 研究工作流——不是 Demo，不是论文附录代码，而是**每天在跑的
生产级工具**。

### 设计哲学

1. **可复现** — 克隆即可运行。打包所有参考数据，零外部数据库依赖。
2. **可分发** — 纯文本 NDJSON 存储，Git 版本控制友好。任何文本编辑器都能打开。
3. **可协作** — 清晰的子项目边界，独立的许可证和文档。每个子项目可以单独使用。
4. **AI 原生** — 工具设计时就考虑了 AI Agent 的操作模式，不是事后补的 API。
5. **经过验证** — 所有工具链都在真实文献分析任务中运行过数百小时，陷阱和边界
   条件已经写进文档。

---

## 推荐技术栈

> ⚡ **经过实操验证的最佳搭配**

### Hermes Agent + DeepSeek V4 Pro

经过数百篇论文的深度分析、数千条知识图谱边的生成、数十次自动化监控任务的运行，
以下组合已被验证为最优解：

| 组件 | 推荐选择 | 官网 |
|------|---------|------|
| 🤖 **AI Agent 宿主** | [**Hermes Agent**](https://github.com/NousResearch/hermes-agent) | [hermes-agent.nousresearch.com](https://hermes-agent.nousresearch.com/docs) |
| 🧠 **大语言模型** | [**DeepSeek V4 Pro**](https://platform.deepseek.com/) | [deepseek.com](https://www.deepseek.com/) |

### 为什么选择这个组合

#### Hermes Agent

Hermes Agent 是由 Nous Research 开发的 CLI AI Agent 框架。与 ChatBot 式的
AI 助手不同，它是一个**工具驱动**的自主代理——能直接操作终端、读写文件、
搜索网页、管理定时任务。

**对 OpenBioMed Lab 的关键价值：**

- 📁 **原生文件系统访问**：直接读写 papers.db、concepts.db、edges.db 等
  NDJSON 知识图谱文件，不需要额外的 API 层
- 🔧 **30+ 内置工具**：`terminal`（执行脚本）、`write_file`（持久化写入）、
  `web_search`（检索）、`cronjob`（定时任务）——全部可在文献工作流中组合使用
- 📚 **Skill 系统**：每个子项目附带 `SKILL.md`，Agent 加载后即掌握完整操作流程，
  无需每次手动解释上下文
- 🕐 **Cron 自动化**：内置定时任务调度器，实现每日文献检索→深度分析→知识图谱
  更新→日报生成的无人值守流水线
- 🔌 **MCP 原生支持**：直接接入 PubMed MCP、arXiv MCP、Fetch MCP 等外部工具，
  与内置 CLI 脚本形成双重保障

> 📖 官网文档：[https://hermes-agent.nousresearch.com/docs](https://hermes-agent.nousresearch.com/docs)
>
> 📦 GitHub：[https://github.com/NousResearch/hermes-agent](https://github.com/NousResearch/hermes-agent)

#### DeepSeek V4 Pro

DeepSeek V4 Pro 是目前性价比最高的推理模型之一。在文献深度分析任务中的表现：

- 📖 **长文本理解**：128K 上下文窗口，可一次加载整篇论文的全文（包括补充材料）
  进行完整 7 层分析
- 🧬 **生物医学知识**：在分子机制推理、通路分析、基因功能解释等任务上表现出色，
  能准确识别修饰位点并绘制因果级联
- 💰 **成本可控**：相比其他同级别模型，API 价格显著更低。一篇完整 S 级 7 层
  分析（T1-T7，含 5 条 C-E-S 链和完整机制级联）的 token 消耗在可控范围内
- 🌐 **中文优势**：对中文生物医学术语的理解准确，适合中英双语文献混合分析

> 📖 官网：[https://www.deepseek.com/](https://www.deepseek.com/)
>
> 🔑 API 平台：[https://platform.deepseek.com/](https://platform.deepseek.com/)

### 实操验证数据

以下数据来自实际运行环境（Hermes Agent + DeepSeek V4 Pro）：

| 指标 | 数值 | 说明 |
|------|------|------|
| 累计分析论文 | 500+ 篇 | S 级 459 篇，知识图谱持续增长 |
| 关联边数量 | 900+ 条 | 5 策略自动生成，纯生物学语义 |
| 每日自动化 | 9:00 HKT | 定时检索 → 分析 → 入库 → 日报 |
| 平均单篇分析 | 完整 7 层 | T1-T7 全层覆盖，无偷懒 |
| 质量自检 | 10 维度 | 每篇写后即时验证空壳检测 |

---

## 子项目

### 1. Literature Learning Suite `v1.3.0`

从研究问题到知识图谱的完整文献认知操作系统。不是文献管理工具——是文献**思考**工具。

[📖 完整文档](literature-learning-suite/README.md) |
[📖 中文指南](literature-learning-suite/GUIDE_ZH.md) |
[📖 English Guide](literature-learning-suite/GUIDE_EN.md)

#### 核心能力

| 能力 | 说明 |
|------|------|
| 📖 **S 级 7 层解剖** | T1 文献档案 → T2 科学问题 → T3 主张-证据链 → T4 机制级联 → T5 隐藏轴 → T6 概念创新 → T7 跨文献关联 |
| 🔗 **5 策略语义关联** | 基于 90,125 基因 + 25,939 通路 + 倒排索引的自动关联边生成 |
| 🏷️ **自动 IF 标注** | 21,800 种期刊的 JCR 2024 影响因子和分区数据，`enrich_paper_if()` 即调即用 |
| 📊 **10 维度质量自检** | 空壳 S 检测、非法边过滤、NDJSON 完整性、全文缓存审计 |
| 🤖 **MCP 集成** | PubMed/arXiv/Fetch/Playwright 四种 MCP 服务器模板和完整设置指南 |
| 🌐 **7 语言文档** | 中文（主）、English、Deutsch、日本語、한국어、Español、العربية |

#### 快速开始

```bash
cd literature-learning-suite
pip install -r scripts/requirements.txt
python scripts/init_workspace.py
```

无需额外配置。打包数据开箱即用。不需要数据库引擎、不需要 API 密钥（PubMed 直连）、不需要 Bioconductor。

#### 一键工作流（配合 Hermes Agent）

在 Hermes Agent 中加载 `literature-learning-suite` skill 后，只需一句话：

```
"检索最近 3 天 Nature/Science/Cell 上与空间转录组学相关的论文，
做 S 级分析，入库，生成日报"
```

Agent 会自动执行：检索 → 去重 → IF 标注 → 全文获取 → 7 层分析 → 写入 papers.db → gen_edges → build_network → gen_digest → selfcheck。

---

## 安装与使用

### 前提条件

```bash
# Python 3.10+
python --version

# Git
git --version

# （可选）Hermes Agent
# 安装见 https://hermes-agent.nousresearch.com/docs
```

### 克隆仓库

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab
```

### 使用子项目

每个子项目是自包含的——进入对应的目录，按照其 README 操作即可：

```bash
cd literature-learning-suite
pip install -r scripts/requirements.txt
python scripts/init_workspace.py
```

### 在 Hermes Agent 中使用

1. 将 `literature-learning-suite/` 复制到 Hermes 的 skills 目录
2. Agent 会自动加载 `SKILL.md` 作为操作协议
3. 配置 MCP 服务器（可选，参见 `assets/templates/mcp-servers.yaml`）
4. 创建定时任务（可选，参见 `references/hermes-monitoring-template.md`）

---

## 许可证

**CC BY-NC-SA 4.0** — 署名-非商业性使用-相同方式共享

| 允许 | 禁止 |
|------|------|
| ✅ 学术研究使用 | ❌ 商业用途 |
| ✅ 个人学习使用 | ❌ 闭源分发衍生作品 |
| ✅ 修改和再分发（需同样许可证） | ❌ 以 SaaS 形式出售 |
| ✅ 大学/研究所内部使用 | ❌ 嵌入商业软件 |

基于本项目二次创作的产品**必须以 CC BY-NC-SA 4.0 或兼容许可证开源**。

📖 完整法律文本：[LICENSE](literature-learning-suite/LICENSE)

---

## 路线图

- [x] Literature Learning Suite v1.3.0 — 文献认知操作系统
- [ ] 病理 AI 分析工具集 — 数字病理图像处理与特征提取
- [ ] 空间转录组学分析流水线 — Xenium/Visium 数据处理
- [ ] 单细胞 RNA-seq 知识图谱 — 细胞类型注释与跨研究整合
- [ ] 多组学数据融合框架 — 基因组+转录组+蛋白质组联合分析

> 欢迎贡献。Fork → 创建子项目文件夹 → 遵循相同结构 → PR。

---

> **OpenBioMed Lab** · CC BY-NC-SA 4.0
>
> Powered by [Hermes Agent](https://github.com/NousResearch/hermes-agent) + [DeepSeek V4 Pro](https://www.deepseek.com/)
