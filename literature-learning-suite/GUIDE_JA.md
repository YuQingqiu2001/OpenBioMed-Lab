# Literature Learning Suite — 完全ユーザーガイド

> **バージョン**: 1.3.0 | **ライセンス**: CC BY-NC-SA 4.0

---

## 目次

1. [概要と設計哲学](#1-概要と設計哲学)
2. [システムアーキテクチャ](#2-システムアーキテクチャ)
3. [クイックスタート](#3-クイックスタート)
4. [コアメソドロジー：S-tier 7層プロトコル](#4-コアメソドロジーs-tier-7層プロトコル)
5. [知識グラフデータモデル](#5-知識グラフデータモデル)
6. [エッジ生成エンジン v3.1](#6-エッジ生成エンジン-v31)
7. [ツールチェーンリファレンス](#7-ツールチェーンリファレンス)
8. [自動モニタリング](#8-自動モニタリング)
9. [トラップライブラリとトラブルシューティング](#9-トラップライブラリとトラブルシューティング)
10. [プラットフォーム別注意事項](#10-プラットフォーム別注意事項)
11. [ベストプラクティス](#11-ベストプラクティス)
12. [FAQ](#12-faq)
13. [用語集](#13-用語集)

---

## 1. 概要と設計哲学

### 1.1 LLSとは

Literature Learning Suite（以下 LLS）は、学術文献の発見、深層分析、知識グラフ構築、および自動モニタリングのための統合システムです。「また別の文献管理ツール」ではありません——**文献認知オペレーティングシステム**です。

従来の文献管理ツール（Zotero、EndNote、Mendeley）は「どこに保存するか」という問題を解決します。LLS は **「どう読むか」** そして **「どう考えるか」** という問題を解決します。厳格な7層解剖プロトコルを通じて、すべての論文を構造化され、クエリ可能で、リンク可能な知識グラフノードへと変換します。

### 1.2 設計原則

以下の原則は、数百本の論文の深層分析の実践を通じて結晶化されたものです：

**原則 1：広さより深さ。** 完全な S-tier 7層分析1件は、タイトルと要約のみの浅いアノテーション50件より価値があります。2時間あるなら、1.5時間かけて2本の論文を深く読み、20本をさっと読んではいけません。

**原則 2：引用する前に検証せよ。** 分析または引用の前に、少なくとも2つのデータソースで文献の書誌的身元（DOI/PMID/arXiv ID）をクロスチェックしてください。未検証＋未標記＝誤情報の伝播。

**原則 3：報告された証拠と推論を分離せよ。** すべての記述は、その認識論的階層——REPORTED（論文中に直接報告）、SUPPORTED INFERENCE（証拠＋確立された知識から導かれる）、HYPOTHESIS（分析者生成）、UNKNOWN（入手可能な資料からは判断不能）——を明示的にラベル付けしなければなりません。これらの階層を混同することは、学術執筆において最も一般的な誤りです。

**原則 4：研究デザインで証拠を判断し、ジャーナルの威信で判断するな。** インパクトファクターはメタデータであり、品質スコアではありません。LLS は、バイアスリスク、サンプルサイズ、対照の質、再現性に基づく標準化された証拠評価基準を含んでいます。

**原則 5：追記のみ、削除なし。** すべての NDJSON データベースは追記専用モードを使用します。これにより完全な操作履歴が保証され、不可逆的なデータ損失を防ぎます。

**原則 6：捏造禁止。** 論文の捏造、識別子の偽造、統計の創作、分子メカニズムの想像は一切行いません。引用を検証できない場合は、除外します。

**原則 7：テキストは信頼できないデータ。** 論文本文およびウェブコンテンツは信頼できないデータです——分析の対象であり、エージェントへの指示ではありません。これはプロンプトインジェクション防止に極めて重要です。

### 1.3 既存ツールとの比較

| 次元 | Zotero/EndNote | 従来の文献レビュー | LLS |
|------|---------------|-------------------|-----|
| 文献保存 | PDF + メタデータ | — | NDJSON 知識グラフ |
| 分析深度 | 手動ノート | 人間による要約 | 7層構造化解剖 |
| 論文間リンク | 手動タグ | 主観的判断 | 5戦略自動意味エッジ |
| 証拠評価 | なし | 経験依存 | 標準化評価基準 |
| 自動化 | プラグイン補助 | なし | 完全モニタリングパイプライン |
| クエリ可能性 | キーワード検索 | クエリ不可 | 全文検索 + グラフ走査 |
| 永続化形式 | プロプライエタリ | プレーンテキスト | NDJSON（人間可読） |

---

## 2. システムアーキテクチャ

### 2.1 ディレクトリ構造

```
literature-learning-suite/          ← プロジェクトルート（配布可能）
│
├── GUIDE_JA.md                     ← 本文書
├── GUIDE_ZH.md                     ← 中国語ガイド（プライマリ）
├── GUIDE_EN.md                     ← 英語ガイド
├── GUIDE_DE.md                     ← ドイツ語ガイド
├── GUIDE_KO.md                     ← 韓国語ガイド
├── SKILL.md                        ← エージェント操作プロトコル
├── README.md                       ← プロジェクト概要
├── LICENSE                         ← CC BY-NC-SA 4.0
├── THIRD_PARTY_DATA.md             ← サードパーティデータ帰属
│
├── scripts/                        ← コアツールチェーン（23+ Python + Node.js + R）
│   ├── init_workspace.py           ← ワークスペース初期化とデータシード
│   ├── literature_search.py        ← マルチソース検索（PubMed/arXiv/Crossref/bioRxiv）
│   ├── search_arxiv.py             ← arXiv 専用検索
│   ├── download_biorxiv_api.py     ← bioRxiv/medRxiv API バッチダウンロード
│   ├── verify_citation.py          ← 文献身元クロス検証
│   ├── normalize_records.py        ← 検索結果の標準化と重複除去
│   ├── validate_records.py         ← JSON Schema 検証
│   ├── fulltext_fetch.py           ← 全文統合ダウンローダー
│   ├── extract_pymupdf.py          ← PDF テキスト/テーブル抽出
│   ├── extract_marker.py           ← PDF OCR 抽出
│   ├── extract_biorxiv_cdp.mjs     ← Chrome CDP 全文抽出（Node.js）
│   ├── biorxiv_chrome_cdp_launcher.bat ← Chrome CDP ランチャー（Windows）
│   ├── kg.py                       ← KG CLI：追加/統計/検索/監査
│   ├── kg_core.py                  ← KG コアライブラリ：CRUD + IF 自動アノテーション
│   ├── ll_common.py                ← 共有ユーティリティ：NDJSON 読み書き、DOI 正規化
│   ├── workspace_paths.py          ← ランタイムパス解決
│   ├── gen_edges.py                ← 意味エッジ生成器 v3.1
│   ├── gen_digest.py               ← 日次ダイジェスト生成
│   ├── build_network.py            ← インタラクティブ力学指向ネットワーク図
│   ├── selfcheck_knowledge_graph.py ← 10次元品質自己チェック
│   ├── export_citations.py         ← BibTeX/CSL エクスポート
│   ├── export_bioc_genes.R         ← 遺伝子/パスウェイ辞書再生器
│   ├── journal_metrics.py          ← ジャーナル指標インポート
│   ├── monitor.py                  ← 再開可能バッチモニター
│   ├── check_assets.py             ← パッケージアセット検証
│   └── requirements.txt            ← Python 依存関係
│
├── assets/
│   ├── data/                       ← 検証済み参照データ（バンドル配布）
│   │   ├── bioc_genes.json         ← 90,125 ヒト＋マウス遺伝子シンボル
│   │   ├── kegg_pathways.json      ← 25,939 KEGG パスウェイ + GO BP 用語
│   │   ├── journal_metrics_2024.json ← 21,800 ジャーナル IF/四分位
│   │   ├── data-manifest.json      ← SHA-256 チェックサム + 出所記録
│   │   ├── evidence-rubric.json    ← 証拠評価基準
│   │   ├── relation-ontology.json   ← エッジタイプオントロジー
│   │   ├── study-designs.json      ← 研究デザイン分類
│   │   ├── search-query-packs.json ← 事前設定検索テンプレート
│   │   └── arxiv-categories.json   ← arXiv カテゴリ体系
│   ├── schemas/                    ← JSON Schema（7種のレコードタイプ）
│   └── templates/                  ← レコードテンプレートと設定テンプレート
│
├── references/                     ← 方法論と操作プロトコル（25+ 文書）
│   ├── data-model.md               ← データモデル仕様
│   ├── deep-analysis-protocol.md   ← 7層分析詳細メソッド
│   ├── s-tier-audit.md             ← 空殻 S 検出ルール
│   ├── s-tier-examples.md          ← 各タイプ論文の分析実例
│   ├── s-tier-upgrade-workflow.md  ← 履歴レコード一括アップグレード手順
│   ├── llm-deep-reasoning-examples.md ← LLM 推論パターンリファレンス
│   ├── gen-edges-v3.md             ← エッジアルゴリズム詳細
│   ├── edge-generation.md          ← エッジタイプリファレンス
│   ├── bioconductor-entity-matching.md ← 遺伝子/パスウェイマッチングロジック
│   ├── full-text-access.md         ← 合法的全文取得ガイド
│   ├── preprint-fulltext.md        ← bioRxiv/medRxiv 抽出
│   ├── self-review-checklist.md    ← コード修正後自己チェックリスト
│   ├── connectivity.md             ← ネットワーク/プロキシ設定
│   ├── cron-troubleshooting.md     ← 無人実行診断
│   ├── automation.md               ← モニタリング自動化設定
│   ├── hermes-monitoring-template.md ← エージェント cron テンプレート
│   ├── mcp-integration.md          ← MCP 統合ガイド
│   ├── mcp-and-tool-routing.md     ← ツールフォールバック戦略
│   ├── journal-metrics-2024.md     ← JCR 指標使用説明
│   ├── pdf-and-ocr.md              ← PDF 抽出方法比較
│   └── bioinfo-tools.md            ← バイオインフォマティクス補助ツール
│
└── tests/                          ← ユニットテスト + 合成データ
```

### 2.2 データフローパイプライン

```
                    研究課題
                       │
                       ▼
              ┌─────────────────┐
              │ 1. 検索戦略立案   │  ← PICO/PECO フレームワーク
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 2. マルチソース検索│  ← PubMed/arXiv/bioRxiv/Crossref
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 3. 身元検証      │  ← デュアルソース DOI/PMID クロスチェック
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 4. 全文取得      │  ← PMC XML / arXiv HTML / CDP 抽出
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 5. S-tier 7層解剖│  ← LLM 深層推論（T1-T7）
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 6. 証拠評価      │  ← バイアスリスク/サンプルサイズ/再現性
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ 7. 知識グラフ永続化│ ← NDJSON 追記
              │  ┌───────────┐  │
              │  │ papers.db │  │
              │  │concepts.db│  │
              │  │ edges.db  │  │
              │  └───────────┘  │
              └────────┬────────┘
                       │
         ┌─────────────┼─────────────┐
         ▼             ▼             ▼
   ┌──────────┐ ┌──────────┐ ┌──────────┐
   │ エッジ生成│ │日次ダイジェスト│ │ ネットワーク│
   │gen_edges │ │gen_digest│ │build_net │
   └──────────┘ └──────────┘ └──────────┘
         │             │             │
         ▼             ▼             ▼
   ┌──────────────────────────────────┐
   │     9. 品質自己チェック + モニタリング    │
   └──────────────────────────────────┘
```

### 2.3 ランタイムワークスペース

`init_workspace.py` によって作成される自己完結型のディレクトリ：

```
my-workspace/                       ← gitignored、ランタイム生成
├── workspace.json                 ← ワークスペースメタデータ
├── papers.db                      ← NDJSON、1行1論文
├── concepts.db                    ← NDJSON、1行1コンセプトノード
├── edges.db                       ← NDJSON、1行1意味エッジ
├── queries.db                     ← NDJSON、検索記録
├── journal_metrics.db             ← NDJSON（オプション、ユーザー定義 IF データ）
├── data/                          ← シードされた遺伝子/パスウェイ辞書のコピー
├── fulltext/                      ← ダウンロードされた全文ドキュメント
├── fulltext_cache/                ← 全文キャッシュ
├── reports/                       ← 生成されたレポート
├── daily_digest/                  ← 日次ダイジェスト
├── config/                        ← モニタリング設定
├── exports/                       ← BibTeX 等のエクスポート
├── imports/                       ← インポート待ちデータ
├── cache/                         ← 汎用キャッシュ
├── biorxiv_api/                   ← bioRxiv API レスポンスキャッシュ
└── logs/                          ← 実行ログ
```

---

## 3. クイックスタート

### 3.1 要件

- Python 3.10+
- pip
- Node.js v24+（bioRxiv CDP 抽出のみ必要）
- R + Bioconductor（遺伝子辞書の再生成時のみ必要）

### 3.2 インストール

```bash
# リポジトリのクローン
git clone <repo-url>
cd literature-learning-suite

# Python 依存関係のインストール
pip install -r scripts/requirements.txt
```

### 3.3 ワークスペースの初期化

```bash
python scripts/init_workspace.py --root ./my-workspace
```

このコマンドは以下を実行します：
1. ワークスペースディレクトリ構造の作成
2. 遺伝子辞書（90,125 シンボル）の `my-workspace/data/` へのコピー
3. パスウェイ/GO 用語（25,939 ラベル）の `my-workspace/data/` へのコピー
4. 空の NDJSON データベースファイル（papers/concepts/edges/queries/journal_metrics.db）の作成
5. デフォルトモニタリング設定のテンプレートからのインストール

### 3.4 環境変数の設定

スクリプトが自動的にワークスペースを見つけられるように：

```bash
# Linux/macOS
export LITERATURE_KG_ROOT="$(pwd)/my-workspace"

# Windows PowerShell
$env:LITERATURE_KG_ROOT = (Resolve-Path ./my-workspace)
```

未設定の場合、スクリプトはデフォルトで `./literature-workspace` を使用します。

### 3.5 最初の検索

```bash
# PubMed 検索
python scripts/literature_search.py pubmed "spatial transcriptomics cancer" -n 10

# arXiv 検索
python scripts/literature_search.py arxiv "single cell foundation model" -n 10

# Crossref 検索
python scripts/literature_search.py crossref "tumor microenvironment review" -n 10

# bioRxiv 検索
python scripts/literature_search.py biorxiv "spatial transcriptomics" -n 10
```

### 3.6 論文の取り込み

```bash
# 検索結果の標準化
python scripts/normalize_records.py search-results.jsonl -o normalized.jsonl

# レコード形式の検証
python scripts/validate_records.py normalized.jsonl

# 知識グラフに追加（自動重複除去）
python scripts/kg.py --root ./my-workspace add normalized.jsonl

# 統計の表示
python scripts/kg.py --root ./my-workspace stats
```

---

## 4. コアメソドロジー：S-tier 7層プロトコル

これが LLS の核心です。すべての研究論文（ジャーナル、IF、プレプリントを問わず）は完全な7層分析を受けなければなりません。

### 4.1 なぜ7層なのか

従来の文献購読は「要約を読む → タグ付け → 2文の要約を書く」で止まります。この浅い処理では：
- 論文の深い論理構造を捕捉できない
- 論文間の暗黙的つながりを発見できない
- 証拠の真の強度を評価できない
- クエリ可能な構造化知識を形成できない

7層プロトコルは、読解プロセスを7つの独立しながらも相補的な認知次元に分解し、分析者（人間またはLLM）に各次元での実質的なアウトプットを強制します。

### 4.2 第1層（T1）：書誌・研究プロファイル

**目標**：論文の検証可能な身元と文脈を確立する。

**内容**：
- 完全なタイトル
- 筆頭著者、責任著者および所属機関
- ジャーナル名（JCR 2024 IF と Q 四分位を自動アノテーション）
- 安定識別子（PMID / DOI / arXiv ID）
- 文献タイプ（Article / Review / Preprint / Clinical Trial / その他）
- 研究デザイン（RCT / コホート / ケースコントロール / 横断 / その他）
- サンプル/モデルシステム、サンプルサイズ
- 全文取得状態（fulltext / abstract_only / metadata_only / unavailable）
- 全文ソース、取得日、バージョン

**自動化度**：この層はツールにより自動抽出可能。`kg_core.enrich_paper_if()` が IF を自動アノテーションします。

### 4.3 第2層（T2）：中核的科学課題

**目標**：論文を1つの反証可能な中核的質問に蒸留し、≥5の検証可能なサブ質問に分解する。

**要件**：
- 中核的質問は**メカニズム的**でなければならない。記述的ではない。「何を見つけたか」ではなく「XがYを通じてZをどう引き起こすか」。
- 各サブ質問は独立して検証可能でなければならない。

**良い例**：
> 中核質問：GDF15は異物受容体シグナルを通じて腫瘍微小環境におけるNK細胞機能不全をどのように誘導するか？
>
> サブ質問：
> 1. どの腫瘍タイプでGDF15が高発現するか？発現レベルと予後の関係は？
> 2. どの受容体がGDF15シグナルを媒介するか？非古典的受容体は関与するか？
> 3. GDF15下流シグナルカスケードはNK細胞の殺傷機能をどう変化させるか？
> 4. この機能不全は可逆的か？薬理学的介入の標的はどこか？
> 5. 生体内モデルでGDF15シグナル遮断は抗腫瘍免疫を回復させるか？

**悪い例（空殻 S）**：
> 中核質問：この論文はGDF15を研究した。
> サブ質問：（空）

### 4.4 第3層（T3）：主張-証拠-統合チェーン（≥5本）

**目標**：論文の中核的発見を ≥5本の独立した C-E-S チェーンに変換し、各チェーンに具体的証拠を含める。

**各チェーンの構造**：

| フィールド | 要件 | 文字数要件 |
|-----------|------|-----------|
| `claim` | 反証可能な結論、自分の言葉で記述 | — |
| `evidence` | 具体的データ：N、効果量、p値、アッセイ系 | **>20文字**（空殻検出閾値） |
| `synthesis` | 証拠が主張をどう支持/弱化するか？分野横断的つながり？ | — |
| `strength` | 1-5★、証拠品質に基づく | — |
| `uncertain` | 代替説明、欠落検証、交絡因子 | — |

**良い例**：
> Claim：PD-L1 1-49%サブグループは術前化学免疫療法から利益を得る。
> Evidence：2847例 NSCLC、SHAP交互作用分析、EFS HR=0.52（95% CI 0.34-0.78）、p=0.002、独立画像評価。
> Synthesis：KEYNOTE-671サブグループ分析が除外したPD-L1「グレーゾーン」集団のエビデンスギャップを埋め、CheckMate-816レジメンがより広い集団に適用可能であることを示唆。
> Strength：★★★★（多施設RCT、大サンプル、独立評価）
> Uncertain：PD-L1検出プラットフォーム間の一致性（22C3 vs 28-8）、アジア人比率の高さが一般化可能性に影響する可能性。

**悪い例（空殻 S）**：
> Claim：新規バイオマーカーを発見した。
> Evidence：著者らは仮説を証明した。
> （evidence フィールド ≤20文字 → 空殻検出トリガー）

### 4.5 第4層（T4）：分子メカニズムカスケード

**目標**：トリガー信号から細胞表現型までの完全な因果連鎖を、修飾部位レベルまで正確にマッピングする。

**構造**：
```
トリガー信号
  │
  ▼
[上流受容体] → [セカンドメッセンジャー/キナーゼ] → [転写因子] → [標的遺伝子] → [細胞表現型]
```

**必須要素**：
- **≥3ステップの因果連鎖**（方向性付き）
- **主要修飾部位**：アミノ酸残基レベルで正確に（例：「NF-κB p65 Ser536 リン酸化」であって「NF-κB活性化」ではない）
- **下流効果**：細胞行動、代謝、相互作用の変化
- **フィードバックループ**：≥1の正または負のフィードバック

**証拠状態ラベル**（各ステップに必須）：
- `demonstrated`：本論文で直接検証
- `supported`：間接的証拠あり
- `background`：確立された背景知識
- `hypothesis`：分析者の推測

**良い例**：
```
トリガー: 腫瘍由来 GDF15
  │
  ▼
GFRAL受容体結合 → [demonstrated: Co-IP, Fig.2A]
  │
  ▼
JAK2-STAT3経路活性化 → [supported: リン酸化抗体, Fig.3B]
  │
  ▼
STAT3 pTyr705 リン酸化 + 核移行 → [demonstrated: 細胞分画+WB, Fig.3C]
  │
  ▼
SOCS3 転写上昇（負のフィードバック） → [demonstrated: qPCR+プロモータールシフェラーゼ]
  │
  ▼
NK細胞殺傷機能低下（CD107a↓, IFN-γ↓） → [demonstrated: フローサイトメトリー, Fig.4]
```

**禁止事項**：
- ❌ 「シグナル経路を活性化」（曖昧すぎる）
- ❌ 修飾部位の捏造
- ❌ 証拠状態の未標記

### 4.6 第5層（T5）：隠れた組織化軸（≥3組）

**目標**：論文が**明示的に述べていない**深層パターン——暗黙の前提、空間的/時間的組織化、選択効果、または実験を貫くロジックを発見する。

**各組の構成**：
- `observation`：論文中の検証可能な具体的事実
- `interpretation`：あなたが発見した深層パターンまたは仮定

**この層はあなたの統合能力をテストするものであり**、論文を書き写す能力ではありません。T5 の observation は論文原文から来なければなりませんが、interpretation は分析者自身が発見しなければなりません。

**良い例**：
> Observation 1：すべての次元削減分析において、腫瘍辺縁サンプルは一貫して腫瘍中心サンプルと別クラスターを形成した。
> Interpretation 1：本研究は暗黙的に「中心」ではなく「辺縁」を疾患決定区画と定義している。これが中心遺伝子シグネチャの予後価値がむしろ低い理由を説明する——中心の分子不均一性が辺縁の微小環境シグナルによって覆い隠されている。
>
> Observation 2：単一細胞分析において、免疫サブセットの変化が腫瘍サブセットの変化に先行して出現した。
> Interpretation 2：研究の時系列データは「免疫優先」の疾患進行モデルを示唆する——微小環境リモデリングが腫瘍進化の原動力であり結果ではない。この解釈が正しければ、早期介入標的は腫瘍細胞ではなく免疫細胞にあるべきである。
>
> Observation 3：薬剤応答者と非応答者の差次的遺伝子において、代謝経路の濃縮度（45%）が免疫経路（12%）を大きく上回った。
> Interpretation 3：論文のタイトルと考察は免疫メカニズムに焦点を当てているが、データの内在的構造は代謝リプログラミングがより根本的な決定因子である可能性を示している。論文のナラティブフレーム（免疫→効果）とデータの重み（代謝→効果）の間に系統的バイアスが存在する。

**禁止事項**：
- ❌ T5 = 論文 Discussion セクションの再掲
- ❌ T5 = 「著者はXを発見した」（これはT3でありT5ではない）

### 4.7 第6層（T6）：概念革新の全体像

**目標**：論文の実質的学術貢献を特定する。要約的な「YにおけるXの役割を研究した」ではない。

**4要素**：

1. **新概念**（≥1）：操作的に定義可能で、独立して引用可能な新概念。漠然とした「新メカニズムを明らかにした」ではない。
2. **覆された/修正された旧説**（≥1）：論文が挑戦または限定した既存の信念は何か？
3. **方法論的ブレイクスルー**（≥1）：論文が貢献した、他者が使用できる新技術能力は何か？
4. **境界条件**：どのような条件下でこの貢献は**成立しない**のか？

**良い例**：
> 新概念：
> - 「イムノメタボリックチェックポイント（Immunometabolic Checkpoint）」：代謝物が免疫受容体を介して媒介する免疫抑制軸と定義され、古典的タンパク質リガンド-受容体免疫チェックポイントとは区別される。
>   - 操作的定義：(a) 小分子代謝物がリガンド、(b) 免疫受容体シグナルを介する、(c) 可逆的免疫細胞機能抑制を引き起こす、の3条件を同時に満たすこと。
>
> 覆された旧説：
> - 「GDF15は単なる食欲抑制因子である」という単一機能仮説を修正し、食欲調節とは独立した腫瘍免疫における役割を証明した。
>
> 方法論的ブレイクスルー：
> - 抗体非依存的な代謝物-受容体結合検出法を確立（SILAC標識 + 化学架橋 + 質量分析）。他のオーファン代謝物受容体の同定に一般化可能。
>
> 境界条件：
> - メカニズムはGDF15高発現腫瘍でのみ成立（中央値の2倍以上）。低発現腫瘍ではこの軸は不活性。
> - 免疫不全モデルで未検証。適応免疫の協調作用を排除できない。

### 4.8 第7層（T7）：論文間関連（≥5本）

**目標**：本論文を他の中核的検証済み文献レコードと、実質的な生物学的関係を用いて接続する。

**各関連の要件**：
- `ref_id`：対象論文の安定識別子（PMID:xxxxx または DOI:10.xxxx/xxxxx）
- `relation`：関係タイプ（下記の有効タイプリスト参照）
- `description`：**60-150語**の説明、関係が存在する理由（WHY）を記述

**有効な関係タイプ**：
- `supports`：本論文の証拠が対象論文の結論を支持
- `contradicts`：本論文の証拠が対象論文と矛盾
- `extends`：本論文が対象論文を有意に拡張
- `replicates`：本論文が対象論文の中核的発見を独立的に再現
- `methodological_complement`：方法論的相補——異なる技術で同一問題を検証
- `shared_mechanism`：分子メカニズムの共有
- `upstream_of` / `downstream_of`：メカニズムの上下流
- `clinical_translation`：基礎から臨床へのトランスレーション
- `shares_disease_model`：同一疾患モデルの使用

**明示的に禁止される非生物学的関係**（`gen_edges.py` 戦略1でフィルター）：
- `same_journal`：同一ジャーナル（生物学的意義なし）
- `same_issue`：同一号
- `same_author`：同一著者
- `same_year`：同一年

### 4.9 空殻 S 検出

S-tier とラベル付けされたレコードが以下の**いずれか**に該当する場合、空殻です：

1. `tier2_subquestions` が空または欠落
2. `tier3_ces_chains` が5本未満
3. いずれかの `evidence` フィールドが ≤20 文字
4. `tier4_mechanism_cascade` のカスケードステップ数が < 3
5. `tier5_hidden_axis` が3組未満
6. `tier7_cross_refs` が5本未満
7. いずれかの T7 relation が禁止リストに該当

**書き込み後検証コマンド**：
```bash
python -c "
import json
with open('./my-workspace/papers.db', 'r', encoding='utf-8-sig') as f:
    for line in f:
        p = json.loads(line)
        if p.get('analysis_tier') != 'S': continue
        chains = p.get('tier3_ces_chains', [])
        empty_ev = sum(1 for c in chains if len(c.get('evidence','')) <= 20)
        print(f'{p[\"id\"]}: T2={len(p.get(\"tier2_subquestions\",[]))} '
              f'T3={len(chains)} empty_evidence={empty_ev} '
              f'T4_steps={len(p.get(\"tier4_mechanism_cascade\",{}).get(\"cascade\",[]))} '
              f'T5={len(p.get(\"tier5_hidden_axis\",[]))} '
              f'T7={len(p.get(\"tier7_cross_refs\",[]))}')
"
```

### 4.10 非研究論文の取り扱い

以下のタイプに該当する論文は S-tier 分析を行いません：
- **ニュース/編集後記** (news/editorial)
- **正誤表/撤回** (erratum/retraction)
- **短報/レター** で実質的な方法/結果内容がないもの

これらの論文には `analysis_tier: "NR"`（Non-Research）とラベル付けし、`core_findings` に実際のタイプ（例：「ニュース記事」「正誤声明」）を記録します。**架空の分析内容の追加は禁止**です。

---

## 5. 知識グラフデータモデル

### 5.1 設計原則

LLS は SQLite や従来のデータベースではなく NDJSON（Newline-Delimited JSON）を使用します。その理由：

1. **人間可読・編集可能**：各行が完全な JSON オブジェクトで、任意のテキストエディタで開ける
2. **バージョン管理フレンドリー**：git diff が行単位で変更を表示
3. **追記専用で不変**：追記専用モードがデータ完全性を保証
4. **CLI でクエリ可能**：`head`、`grep`、`jq` などの標準ツールが直接使用可能
5. **依存関係ゼロ**：データベースエンジンのインストール不要

### 5.2 papers.db

論文マスターデータベース。各行に1件の完全な論文レコード。

**必須フィールド**：`id`（安定識別子）、`title`、`source`（検索ソース）、`retrieved_at`

**S-tier 分析フィールド**（v4.0 標準）：
```json
{
  "id": "PMID:42251595",
  "title": "完全なタイトル",
  "authors": ["筆頭著者", "..."],
  "journal": "ジャーナル正式名",
  "impact_factor": 63.1,
  "journal_quartile": "Q1",
  "doi": "10.xxxx/xxxxx",
  "pmid": "42251595",
  "source": "pubmed",
  "retrieved_at": "2026-06-09T09:00:00",
  "fulltext_status": "fulltext",
  "analysis_tier": "S",
  "analysis_method": "LLM_deep_reasoning_S_tier_v4.0",
  "tier2_core_question": "...",
  "tier2_subquestions": ["Q1", "Q2", "Q3", "Q4", "Q5"],
  "tier3_ces_chains": [{"chain_id":1, "claim":"...", "evidence":"...", "synthesis":"...", "strength":3, "uncertain":"..."}],
  "tier4_mechanism_cascade": {"trigger":"...", "cascade":["Step1","Step2","Step3"], "key_modifications":[...], "downstream_effects":"...", "feedback":[...], "evidential_status":{...}},
  "tier5_hidden_axis": [{"observation":"...", "interpretation":"..."}],
  "tier6_concept_innovation": {"new_concepts":[...], "overturned_views":[...], "methodological_breakthroughs":[...], "boundary_conditions":"..."},
  "tier7_cross_refs": [{"ref_id":"PMID:xxxxx", "relation":"extends", "description":"..."}],
  "entities": {"genes":["GDF15"], "pathways":["JAK-STAT"], "cell_types":["NK cells"], "diseases":["cancer"]}
}
```

### 5.3 concepts.db

コンセプト/エンティティデータベース。

```json
{
  "id": "CONCEPT:immunometabolic_checkpoint",
  "name": "イムノメタボリックチェックポイント",
  "type": "mechanism",
  "definition": "小分子代謝物が免疫受容体を介して媒介する免疫抑制軸",
  "source_papers": ["PMID:42251595", "PMID:39988000"],
  "created_at": "2026-06-09T09:00:00"
}
```

`type` の有効値：`mechanism`, `disease`, `method`, `cell_type`, `pathway`, `drug`, `gene`, `phenomenon`, `hypothesis`

### 5.4 edges.db

意味エッジデータベース。

```json
{
  "source": "PMID:42251595",
  "target": "PMID:39988000",
  "relation": "extends",
  "description": "PMID:42251595は PMID:39988000 が同定した GDF15-GFRAL 結合に基づき、GFRAL 下流の JAK2-STAT3-SOCS3 負のフィードバックループをさらに発見し、GDF15 シグナルをリガンド-受容体認識から完全なシグナル伝達カスケードへと拡張した。両者は共同して GDF15 免疫抑制軸の分子的枠組みを構成する。",
  "provenance": "analyst",
  "created_at": "2026-06-09T10:00:00"
}
```

### 5.5 queries.db

検索記録。監査と再現のために使用。

```json
{
  "source": "pubmed",
  "query": "spatial transcriptomics cancer",
  "executed_at": "2026-06-09T09:00:00",
  "result_count": 42,
  "parameters": {"retmax": 50, "sort": "date"}
}
```

---

## 6. エッジ生成エンジン v3.1

### 6.1 設計目標

`gen_edges.py` の目標は、論文間に**生物学的に意味のある意味的関連**を確立することです。従来のキーワード共起や引用グラフベースのアプローチとは異なり、v3.1 は5つの相補的戦略を使用し、すべて転置インデックスで実装されています（時間計算量 O(M)、O(N²) ではありません）。

### 6.2 5戦略の詳細

#### 戦略 1：明示的参照
- **データソース**：論文の `tier7_cross_refs` フィールド（LLM 深層分析の結果）
- **フィルター**：非生物学的関係（`same_journal`, `same_issue`, `same_author`）を除外
- **特徴**：最高品質（人間/LLM 専門家判断）、ただしカバレッジは分析深度に依存

#### 戦略 2：共有分子（≥2共有）
- **データソース**：90,125 ヒトおよびマウス遺伝子シンボル（Bioconductor エクスポート）
- **マッチングロジック**：二段階マッチング——
  1. 正規表現で論文テキストから候補語を抽出
  2. 遺伝子セット内で O(1) ルックアップ確認
- **利点**：高精度（ML/TNF のような略語の偽陽性を排除）
- **制限**：1遺伝子あたり最大15論文（TNF のような高頻度遺伝子が過剰なエッジを生成するのを防止）

#### 戦略 2.5：テキスト重複（≥4共有キーワード）
- **データソース**：論文のコア発見（`core_findings` または T3 claims）
- **マッチングロジック**：
  1. Bag-of-words を抽出（ストップワード除去）
  2. 単語→論文の転置インデックスを構築
  3. 各論文ペアの共有単語数を集計
- **特徴**：主力戦略。全エッジの約70%に貢献

#### 戦略 3：同一疾患×同一手法
- **データソース**：論文の `diseases` および `methods`/`technologies` フィールド
- **マッチングロジック**：疾患ラベル × 手法ラベルのクロス積
- **特徴**：粒度は粗いが、カバレッジは低く精度は高い

#### 戦略 4：隠れた軸
- **データソース**：論文の `tier5_hidden_axis` フィールド
- **マッチングロジック**：T5 から深層パターンのキーワード（paradigm, bias, survivor, selection など）を抽出し、上位200論文に対して共鳴マッチング
- **特徴**：最も深い関連タイプ。暗黙的パラダイムレベルでの共通性を捕捉

#### 戦略 5：コンセプトノード
- **データソース**：`concepts.db`
- **マッチングロジック**：各コンセプトの `source_papers` に対して `defines_concept` エッジを生成
- **特徴**：コンセプトノードを論文ノードに接続

### 6.3 実行コマンド

```bash
# 重要！毎回実行前にバイトコードキャッシュをクリアすること
rm -rf scripts/__pycache__

# エッジ生成器を実行（-B で新規 .pyc の書き込みを禁止）
python -B scripts/gen_edges.py

# インタラクティブネットワーク図を更新
python scripts/build_network.py
```

### 6.4 出力形式

```json
{
  "source": "PMID:42251595",
  "target": "PMID:39988000",
  "relation": "shares_molecules",
  "description": "共有分子: genes:GDF15, genes:TNF, pathways:JAK-STAT signaling",
  "metadata": {"shared_entities": ["genes:GDF15", "genes:TNF", "pathways:JAK-STAT signaling"]},
  "provenance": "deterministic_molecule_index",
  "created_at": "2026-06-09T10:00:00"
}
```

---

## 7. ツールチェーンリファレンス

### 7.1 検索と発見

#### `literature_search.py`
マルチソース学術文献検索。

```
使用法: literature_search.py <source> <query> [options]

ソース (source):
  pubmed    - PubMed E-utilities API
  arxiv     - arXiv API
  crossref  - Crossref API
  biorxiv   - bioRxiv API

オプション:
  -n N      - 最大返却数（デフォルト 20）
  -y YEAR   - 年フィルター
  -o FILE   - 出力ファイル（デフォルト stdout）
```

#### `search_arxiv.py`
arXiv 専用検索。著者、カテゴリ、日付範囲フィルター対応。

```
使用法: search_arxiv.py [--author NAME] [--category CAT] [--max N]
```

#### `download_biorxiv_api.py`
bioRxiv/medRxiv 日付別バッチダウンロード。

```
使用法: download_biorxiv_api.py --date YYYY-MM-DD [--source biorxiv|medrxiv]
```

### 7.2 検証と標準化

#### `verify_citation.py`
文献身元クロス検証。

```
使用法: verify_citation.py --doi 10.xxxx/xxxxx
       verify_citation.py --pmid 12345678
       verify_citation.py --arxiv 2401.01234
```

検証フロー：
1. 主要データソースでレコードの存在を確認
2. タイトル、著者、年、ジャーナルの一貫性をチェック
3. 撤回/正誤/懸念表明をチェック
4. 重要レコードは第2データソースでクロス検証

#### `normalize_records.py`
検索結果の標準化と重複除去。

重複除去の優先順位：PMID > arXiv ID > 正規化 DOI > 正規化タイトル+年

### 7.3 全文取得

#### `fulltext_fetch.py`
統合全文ダウンローダー。最適なパスを自動選択。

```
パス優先順位：
1. PubMed Central OA XML（PMCID が必要）
2. arXiv HTML（認証不要）
3. bioRxiv/medRxiv API 要約（常に利用可能）
4. ユーザー提供のローカル PDF
```

#### `extract_biorxiv_cdp.mjs`
Chrome DevTools Protocol 全文抽出器。Cloudflare 保護された bioRxiv/medRxiv ページを処理。

```
前提条件：
1. Chrome をリモートデバッグモードで起動（ポート 9223）
2. ブラウザで手動でセキュリティ検証を完了（自動化不可）
3. その後、本スクリプトを実行してレンダリング済み本文を抽出

使用法: node scripts/extract_biorxiv_cdp.mjs --doi 10.1101/XXXX --port 9223
```

**設計原則**：自動 CAPTCHA 解法やアクセス制御の回避は使用しません。CDP 抽出器は「手動選択→コピー→貼り付け」の手動作業を置き換えるだけです。

### 7.4 知識グラフ操作

#### `kg.py`
知識グラフコマンドラインツール。

```
使用法: kg.py --root <workspace> <command> [args]

コマンド:
  add <file>     - 論文レコードを追加（JSONL または NDJSON）、自動重複除去
  stats          - 統計を表示（各DBレコード数、ソース分布、全文状態）
  search <query> - 全文検索（AND 論理、大文字小文字区別なし）
  audit          - 完全性監査（重複 ID、必須フィールド欠落、JSONパースエラー）
```

#### `kg_core.py`
知識グラフコアライブラリ（Python API）。

```python
from kg_core import (
    add_paper,          # 論文追加（自動重複除去+IFアノテーション）
    enrich_paper_if,    # JCR 2024 IF 自動アノテーション
    get_stats,          # 統計取得
    get_recent_papers,  # 最近N日間の論文取得
    search_papers,      # 論文検索
    lookup_journal_impact_factor,  # ジャーナル IF 検索
    write_daily_digest, # 日次ダイジェスト書き込み
)

# 例
paper = {'title': '...', 'pmid': '12345', 'journal': 'Nature', 'source': 'pubmed'}
paper = enrich_paper_if(paper)  # IF=50.5, Q1 を自動補完
added = add_paper(paper)  # True（新規論文）または False（重複）
```

### 7.5 分析と出力

#### `gen_edges.py`
エッジ生成器 v3.1（第6節で詳述）。

#### `gen_digest.py`
日次ダイジェスト生成器。

```
生成内容：
- 論文総数と S/A/B レベル分布
- エッジタイプ分布
- S-tier 論文リスト（上位15件、ジャーナルとコア発見付き）

出力：daily_digest/YYYY-MM-DD.md
```

#### `build_network.py`
インタラクティブ HTML ネットワーク図生成器。

```
出力：network.html（ダブルクリックで開く）

可視化ルール：
- ノード色：PubMed=緑, arXiv=赤, bioRxiv=青, medRxiv=水色, コンセプト=黄
- エッジ色：論文間推論=青, 共有遺伝子=オレンジ, 共有パスウェイ=紫, 論文→コンセプト=黄
- ノードサイズ ∝ claims 数
- 白枠 = 完全深層分析済み
- 完全オフライン、外部依存なし
```

#### `selfcheck_knowledge_graph.py`
10次元品質自己チェック。

| 次元 | チェック内容 |
|------|-------------|
| ファイル一覧 | 全ディレクトリスキャン、ファイル数/サイズ統計 |
| 禁止残留物 | chrome_cdp_profile、cookie ファイル、一時テストファイル |
| CDP ポート | 9222/9223 がまだ開いているか |
| DB 完全性 | NDJSON パースエラー、重複 ID、必須フィールド欠落 |
| S-tier 品質 | 空殻 S / 弱 S 検出（1論文あたり7項目） |
| コンセプト監査 | 重複 ID、name 欠落 |
| エッジ監査 | 不正な非生物関係、説明なしエッジ、孤立エッジ、自己ループ、重複エッジ |
| 全文キャッシュ | 命名規則、過小ファイル、Cloudflare 残留、重複内容 |
| エッジ統計 | 戦略別エッジ数分布 |
| ネットワーク一貫性 | クロスリファレンスの有効性 |

#### `export_citations.py`
BibTeX/CSL 形式へのエクスポート。

```
使用法: export_citations.py <papers.db> --format bibtex -o library.bib
```

---

## 8. 自動モニタリング

### 8.1 モニタリング設定

モニタリングジョブは `config/monitor-job.json` で定義します：

```json
{
  "name": "日次トップジャーナル文献深層学習",
  "schedule": "0 9 * * *",
  "sources": ["pubmed", "arxiv", "biorxiv"],
  "date_window_days": 1,
  "max_papers_per_source": 50,
  "analysis_tier": "S",
  "dedup": true,
  "resumable": true
}
```

### 8.2 モニタリングの実行

```bash
python scripts/monitor.py ./my-workspace/config/monitor-job.json --root ./my-workspace
```

モニターの特性：
- **再開可能**：中断後、前回の続きから再開
- **バッチ処理**：少量ずつ処理し、メモリオーバーフローを防止
- **冪等**：繰り返し実行しても重複取り込みなし（ID ベースの重複除去）
- **障害分離**：単一論文の分析失敗が全体のフローに影響しない

### 8.3 エージェント cron タスク

Hermes Agent 環境で使用する場合、cron タスクを作成できます：

```json
{
  "name": "日次文献深層学習",
  "schedule": "0 9 * * *",
  "skills": ["literature-learning-suite"],
  "enabled_toolsets": ["terminal", "file", "web"],
  "workdir": "/path/to/literature-learning-suite"
}
```

**重要**：cron タスクではインタラクティブツールの使用を禁止してください。bioRxiv CDP 抽出はユーザーの手動操作が必要であり、無人実行には適しません。

### 8.4 タイムゾーン設定

タスクのトリガー時刻が不正確な場合、エージェント設定の `timezone` 設定を確認してください。未設定の場合、スケジューラーは UTC を使用する可能性があります。

---

## 9. トラップライブラリとトラブルシューティング

以下は本番環境で実際に発生したエラーを深刻度順に並べたものです。

### TRAP-001：.pyc バイトコードゴースト
- **現象**：ソースコードを修正しても動作が変わらない
- **根本原因**：Python が `__pycache__/` 内の古い .pyc ファイルをロードしている
- **予防**：`gen_edges.py` の修正後は毎回 `rm -rf scripts/__pycache__ && python -B` を実行
- **発生シーン**：`gen_edges.py` の戦略コード修正時に最も頻発

### TRAP-002：サンドボックス書き込み不可視
- **現象**：`execute_code` で書き込んだファイルがディスク上に存在しない
- **根本原因**：一部のエージェント環境が一時サンドボックスを使用し、呼び出し終了後にファイルを破棄
- **予防**：永続化書き込みには `write_file` ツールまたは `terminal` + Python を使用
- **パターン**：`write_file(path, content)` → `terminal("python path/to/script.py")`

### TRAP-003：.db ファイル拒否
- **現象**：`read_file` が `Cannot read binary file 'papers.db'` エラー
- **根本原因**：`.db` 拡張子がバイナリファイルとして認識されるが、実際のファイルは NDJSON テキスト
- **予防**：常に `terminal` + `python -c` または `head` で .db ファイルを読み取る

### TRAP-004：Unicode 文字による SyntaxError
- **現象**：`SyntaxError: invalid character '→' (U+2192)`
- **根本原因**：矢印（→）、曲線引用符（""）、ギリシャ文字（ΔΨ）が一部のシェルで誤解釈される
- **予防**：Python スクリプト内で非 ASCII 文字を避ける。置換ルール：
  - `→` → `->`
  - `""` → `'`
  - `ΔΨ` → `Delta`/`Psi`

### TRAP-005：arXiv API サイレント失敗
- **現象**：0件の結果が返るが、エラーメッセージなし
- **根本原因**：HTTPS で `export.arxiv.org` にアクセスすると 301 リダイレクトが返り、`-L` なしではサイレント失敗
- **予防**：`curl -sL "http://export.arxiv.org/api/query?..."` を使用（HTTP + -L に注意）

### TRAP-006：PubMed プロキシ経由で遅延
- **現象**：PubMed クエリがタイムアウトまたは極端に遅い
- **根本原因**：PubMed API は一部の地域で直結の方が高速（~780ms）。プロキシ経由だと逆に遅くなる
- **予防**：PubMed クエリ前に `unset http_proxy https_proxy`

### TRAP-007：空殻 S-tier 論文
- **現象**：`analysis_tier: "S"` だが T2-T7 フィールドが空
- **根本原因**：ラベルのみ変更し、実質的内容を書いていない
- **検出**：`selfcheck_knowledge_graph.py` を実行し、`s_empty` カウントを確認
- **予防**：毎回書き込み後に書き込み後検証スクリプトを実行（4.9節参照）

### TRAP-008：非生物エッジ汚染
- **現象**：`gen_edges.py` が `same_journal` エッジを生成
- **根本原因**：LLM が書いた T7 cross_refs に非生物学的関係が含まれる可能性
- **予防**：戦略1に `NON_BIO_RELS` フィルターを内蔵
- **検証**：`selfcheck_knowledge_graph.py` のエッジ監査次元

### TRAP-009：bioRxiv JATS XML 403
- **現象**：`curl https://www.biorxiv.org/content/10.1101/XXXX.source.xml` が 403 を返す
- **根本原因**：bioRxiv/medRxiv がソース XML および PDF エンドポイントへのプログラム的アクセスをブロック（2026年）
- **代替手段**：Chrome CDP 抽出（`extract_biorxiv_cdp.mjs`）または API 要約（300-500語）を使用

### TRAP-010：ジャーナル IF ファジーマッチ誤判定
- **現象**："Cell" の IF が "Cell Reports" の IF とマッチングされる
- **根本原因**：部分文字列マッチングの粒度が粗すぎる
- **修正**：`kg_core.py` に長さ比閾値 0.7 を追加（`ratio = min(len(norm),len(key)) / max(len(norm),len(key))`）

---

## 10. プラットフォーム別注意事項

### Linux / macOS

```bash
# コマンドプレフィックス
python3 scripts/init_workspace.py

# Chrome CDP 起動
google-chrome --remote-debugging-port=9223
# または
chromium --remote-debugging-port=9223

# 遺伝子辞書再生（R + Bioconductor が必要）
Rscript scripts/export_bioc_genes.R ./my-workspace
```

### Windows

```bash
# git-bash (MSYS) または WSL の使用を推奨
# コマンドプレフィックス
python scripts/init_workspace.py   # python3 ではない

# Chrome CDP 起動
./scripts/biorxiv_chrome_cdp_launcher.bat

# パス形式
/c/Users/.../my-workspace   # MSYS スタイル
C:\path\to\my-workspace   # Windows ネイティブスタイル
# 両方とも使用可能
```

**注意**：Windows MSYS bash で `python -c "..."` インラインコードを実行する場合、文字列内での非 ASCII 文字の使用を避けてください（TRAP-004 参照）。

### クロスプラットフォーム共通

- バンドルされた遺伝子/パスウェイ辞書と JCR ジャーナル指標はそのまま使用可能
- `.db` ファイルはプレーンテキスト NDJSON であり、任意のエディタで開ける
- すべてのスクリプトは `pathlib.Path` を使用し、パス区切り文字を自動処理
- `enrich_paper_if()` が21,800ジャーナルのデータベースから IF を自動アノテーション

---

## 11. ベストプラクティス

### 11.1 検索戦略

- **狭く始めて、徐々に広げる**：最初に精密な用語で高精度の結果を得てから、徐々に範囲を拡大
- **PubMed を最初に、arXiv を後に**：PubMed の応答が高速（直結 ~780ms）、カバレッジも広い
- **すべての検索を記録**：完全なクエリ文字列、実行時刻、結果数——再現性と監査に極めて重要
- **ゼロ結果を「存在しない」と解釈しない**：まずクエリ構文、ネットワーク接続、日付フィルター、レート制限を確認

### 11.2 分析戦略

- **全文を優先**：全文あり→完全な7層分析を実行。要約のみ→主張の強度を制限し、`abstract_only` とラベル付け
- **バッチ処理禁止**：論文は1本ずつ独立して分析する。S-tier 7層分析は1本ずつ。前の論文のテンプレートをコピー＆ペーストすることが「空殻 S」の最大の原因
- **証拠は具体的に**：「2847例 NSCLC、SHAP 分析、EFS HR=0.52」であって「著者が証明した」ではない
- **メカニズムは部位レベルで正確に**：「NF-κB p65 Ser536 リン酸化」であって「シグナル経路を活性化」ではない

### 11.3 メンテナンス戦略

- **毎日 `selfcheck_knowledge_graph.py`** を cron タスクの最終ステップとして実行。日報の前に必ず監査
- **定期的に .pyc をクリア**：特に `gen_edges.py` の修正後
- **週次バックアップ**：papers.db、edges.db、concepts.db は最も重要な3ファイル
- **ログの監視**：`logs/` ディレクトリの異常を定期的にチェック

### 11.4 チーム協業

- **1人1ワークスペース**：各研究者が自分のワークスペースを維持し、同時書き込みの競合を回避
- **コンセプトライブラリの共有**：concepts.db はワークスペース間でマージ可能
- **エッジは再構築可能**：edges.db は完全に gen_edges.py が生成するため、いつでも再構築可能
- **Git の使用**：.db ファイル（NDJSON テキスト）はバージョン管理可能。git diff が変更を表示

---

## 12. FAQ

**Q: なぜ .db ファイルは SQLite ではなくテキストなのですか？**
A: NDJSON は人間可読、バージョン管理フレンドリー、CLI クエリ可能、外部依存ゼロという利点があります。文献知識グラフのような「一度書き、何度も読み、時々追記する」ワークロードには、NDJSON が SQLite より適しています。

**Q: 増分更新はどのように処理しますか？**
A: `add_paper()` と `kg.py add` には安定 ID に基づく重複除去ロジックが組み込まれています。新しいレコードを追記するだけで、既存のレコードは上書きされません。

**Q: gen_edges.py はどのくらいの頻度で実行すべきですか？**
A: 新しい論文が追加されるたびに実行します。日次モニタリングでは、cron タスクの最終ステップとして実行します（.pyc クリア → gen_edges → build_network → selfcheck）。

**Q: 独自のジャーナル IF データを使用できますか？**
A: はい。IF データを `my-workspace/journal_metrics.db`（NDJSON 形式、フィールドは assets/data/journal_metrics_2024.json と一致）に配置してください。システムはワークスペースのデータを優先し、次にバンドルデータにフォールバックします。

**Q: bioRxiv の全文をどうやって取得しますか？**
A: API は 300-500 語の拡張要約を提供します（常に利用可能）。JavaScript レンダリングされた全文ページは Chrome CDP 抽出（`extract_biorxiv_cdp.mjs`）が必要で、このパスでは Cloudflare 検証のためにブラウザの手動操作が必要です。

**Q: 日本語文献はサポートされていますか？**
A: 現在対応している検索ソース（PubMed/arXiv/bioRxiv/Crossref）は主に英語文献向けです。日本語文献は手動でレコードを作成し、`normalize_records.py` + `kg.py add` で取り込むことが可能です。分析プロトコル自体は言語に依存しません。

**Q: 古いバージョンの分析レコードをアップグレードするには？**
A: `references/s-tier-upgrade-workflow.md` を参照してください。基本的な流れ：監査 → 一括アップグレード（10件/セット）→ エッジとネットワーク図の再構築 → 再監査。

---

## 13. 用語集

| 日本語 | English | 説明 |
|--------|---------|------|
| 知識グラフ | Knowledge Graph | NDJSON 形式で保存された論文-コンセプト-エッジのネットワーク |
| S-tier 分析 | S-tier Analysis | 完全な7層解剖プロトコル分析 |
| 空殻 S | Empty-shell S | S とラベル付けされているが T2-T7 に実質的内容がないレコード |
| 主張-証拠-統合チェーン | Claim-Evidence-Synthesis Chain (CES) | T3 層のコアユニット |
| 隠れた組織化軸 | Hidden Organizing Axis | T5 層。論文が明示的に述べていない深層パターン |
| 意味エッジ | Semantic Edge | 生物学的に意味のある論文間関係 |
| NDJSON | Newline-Delimited JSON | 各行に1つの完全な JSON オブジェクトを含む保存形式 |
| 転置インデックス | Inverted Index | gen_edges.py のコアデータ構造（O(M) 計算量） |
| ワークスペース | Workspace | ランタイムデータディレクトリ（my-workspace/） |
| 安定識別子 | Stable Identifier | PMID/DOI/arXiv ID など、永続的に引用可能な ID |
| 証拠評価基準 | Evidence Rubric | バイアスリスク/サンプルサイズ/再現性に基づく評価体系 |
| JCR | Journal Citation Reports | Clarivate 社のジャーナル引用レポート |
| CDP | Chrome DevTools Protocol | ブラウザ自動化プロトコル |
| 二段階マッチング | Two-stage Matching | 正規表現候補抽出 → セットルックアップ確認 |
| イムノメタボリックチェックポイント | Immunometabolic Checkpoint | 代謝物が免疫受容体を介して媒介する免疫抑制軸 |
| 認識論的階層 | Epistemological Tier | 記述の確からしさのレベル（REPORTED/SUPPORTED INFERENCE/HYPOTHESIS/UNKNOWN） |
| PICO/PECO フレームワーク | PICO/PECO Framework | 検索戦略立案のための標準化フレームワーク |

---

> 本文書は Literature Learning Suite プロジェクトによりメンテナンスされています。
> バージョン 1.3.0、最終更新 2026-06-09。
> ライセンス：CC BY-NC-SA 4.0（表示-非営利-継承）
