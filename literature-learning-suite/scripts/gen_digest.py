"""
每日文献深度学习日报生成器
=========================
从 kg_core 动态读取数据，生成 markdown 日报。
"""

import json
from datetime import date
from pathlib import Path
from collections import Counter

from workspace_paths import KG_ROOT

PAPERS_DB = KG_ROOT / "papers.db"
EDGES_DB = KG_ROOT / "edges.db"
DIGEST_DIR = KG_ROOT / "daily_digest"


def load_ndjson(path):
    records = []
    if path.exists():
        with open(path, 'r', encoding='utf-8-sig') as f:
            for line in f:
                if not line.strip(): continue
                try: records.append(json.loads(line))
                except: pass
    return records


def analysis_summary(paper):
    analysis = paper.get('analysis')
    analysis = analysis if isinstance(analysis, dict) else {}
    chains = paper.get('tier3_ces_chains') or analysis.get('ces_chains') or []
    core = (
        paper.get('core_findings')
        or analysis.get('core_question')
        or next(
            (
                chain.get('claim', '')
                for chain in chains
                if isinstance(chain, dict) and chain.get('claim')
            ),
            '',
        )
    )
    tier = paper.get('analysis_tier') or ('S' if chains else '?')
    return tier, core


def main():
    today = date.today().isoformat()
    papers = load_ndjson(PAPERS_DB)
    edges = load_ndjson(EDGES_DB)

    # 统计
    tiers = Counter(analysis_summary(p)[0] for p in papers)
    s_papers = [p for p in papers if analysis_summary(p)[0] == 'S']
    sources = Counter(p.get('source', 'unknown') for p in papers)
    edge_types = Counter(e.get('relation', '?') for e in edges)

    digest = f"""# 每日文献深度学习日报 — {today}

## 执行摘要
- **论文总量**: {len(papers)} 篇
- **S级**: {tiers.get('S', 0)} | A级: {tiers.get('A', 0)} | B级: {tiers.get('B', 0)}
- **关联边**: {len(edges)} 条
- 来源: {', '.join(f'{k} {v}' for k,v in sources.most_common())}

## 关联边类型分布
| 类型 | 数量 |
|------|------|
"""
    for rel, cnt in edge_types.most_common():
        digest += f"| {rel} | {cnt} |\n"

    if s_papers:
        digest += "\n## S 级论文 ({}) \n\n".format(len(s_papers))
        for p in s_papers[:15]:
            title = (p.get('title') or 'Untitled')[:80]
            journal = p.get('journal', '?')
            findings = analysis_summary(p)[1][:200]
            digest += f"### {title}\n"
            digest += f"- 期刊: {journal}\n"
            if findings:
                digest += f"- 核心发现: {findings}\n"
            digest += "\n"

    digest += f"\n---\n*由 Hermes Agent gen_digest.py 自动生成*\n"

    out_path = DIGEST_DIR / f"{today}.md"
    DIGEST_DIR.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(digest)

    print(f"日报已生成: {out_path}")


if __name__ == "__main__":
    main()
