"""
文献知识图谱核心模块
=====================
NDJSON-based knowledge graph for biomedical literature.
Papers indexed with entities (genes, pathways, cell types, diseases),
relevance scores, and citation edges.

Storage layout under $LITERATURE_KG_ROOT:
    papers.db         - NDJSON, 每行一篇论文
    concepts.db       - NDJSON, 每行一个概念/实体
    edges.db          - NDJSON, 每行一条关系边
    daily_digest/     - 每日简报 markdown
"""

import json
import hashlib
from datetime import datetime, date, timedelta
from pathlib import Path
import re

from workspace_paths import KG_ROOT, JOURNAL_METRICS_PATH, BUNDLED_JOURNAL_METRICS, ensure_workspace

# 懒加载期刊指标
_journal_metrics: dict = None  # {normalized_name: {jif, quartile, category, ...}}
_journal_metrics_loaded: bool = False
PAPERS_DB = KG_ROOT / "papers.db"
CONCEPTS_DB = KG_ROOT / "concepts.db"
EDGES_DB = KG_ROOT / "edges.db"
DIGEST_DIR = KG_ROOT / "daily_digest"


def _normalize_name(name: str) -> str:
    """标准化期刊名用于匹配"""
    if not name:
        return ""
    return re.sub(r'[^a-z0-9]', '', name.lower().strip())


def _safe_float(val):
    """安全转换 JIF 值，处理 '<0.1' 等特殊情况"""
    if val is None:
        return 0.0
    if isinstance(val, (int, float)):
        return float(val)
    s = str(val).strip()
    if not s:
        return 0.0
    # 处理 '<0.1', '>100' 等情况
    s = re.sub(r'[<>]', '', s)
    try:
        return float(s)
    except ValueError:
        return 0.0


def _load_journal_metrics():
    """懒加载期刊指标数据库，优先用户数据，回退到打包数据"""
    global _journal_metrics, _journal_metrics_loaded
    if _journal_metrics_loaded:
        return
    _journal_metrics = {}
    # Try user path first; if no valid entries, try bundled data
    for source_path in (JOURNAL_METRICS_PATH, BUNDLED_JOURNAL_METRICS):
        if not source_path.exists():
            continue
        text = source_path.read_text(encoding="utf-8-sig").strip()
        if not text:
            continue
        if text.startswith("["):
            data = json.loads(text)
        else:
            data = [
                json.loads(line)
                for line in text.splitlines()
                if line.strip()
            ]
        if not data:
            continue
        _journal_metrics = {}
        for entry in data:
            if not isinstance(entry, dict) or not entry.get("name"):
                continue
            key = _normalize_name(entry["name"])
            _journal_metrics[key] = {
                "name": entry["name"],
                "abbr": entry.get("abbr_name", ""),
                "jif": _safe_float(entry.get("jif", 0)),
                "jif_5y": _safe_float(entry.get("jif_5y", 0)),
                "category": entry.get("category_name", ""),
                "category_detail": entry.get("category_raw", ""),
                "quartile": entry.get("quartile", ""),
                "rank": entry.get("rank", ""),
                "rank_total": entry.get("rank_total", ""),
            }
        if _journal_metrics:
            break  # User data loaded successfully
    _journal_metrics_loaded = True


def lookup_journal_impact_factor(journal_name: str) -> dict:
    """
    查找期刊影响因子。
    返回: {"found": bool, "jif": float, "name": str, "quartile": str, ...}
    """
    _load_journal_metrics()
    
    if not journal_name or not _journal_metrics:
        return {"found": False, "jif": 0, "name": journal_name}
    
    norm = _normalize_name(journal_name)
    
    # 精确匹配
    if norm in _journal_metrics:
        return {"found": True, **_journal_metrics[norm]}
    
    # 模糊匹配 — 仅当部分匹配且长度比 > 0.7（防止 "Cell" 误匹配 "Cell Reports"）
    for key, info in _journal_metrics.items():
        if norm in key or key in norm:
            ratio = min(len(norm), len(key)) / max(len(norm), len(key))
            if ratio > 0.7:
                return {"found": True, **_journal_metrics[key]}
    
    return {"found": False, "jif": 0, "name": journal_name}


def enrich_paper_if(paper: dict) -> dict:
    """为论文补充影响因子信息"""
    journal = paper.get("journal", "")
    if not journal:
        return paper
    
    result = lookup_journal_impact_factor(journal)
    if result["found"]:
        paper["impact_factor"] = result["jif"]
        paper["impact_factor_5y"] = result.get("jif_5y", 0)
        paper["journal_quartile"] = result.get("quartile", "")
        paper["journal_category"] = result.get("category", "")
    
    return paper


def get_top_journals(limit: int = 20) -> list:
    """获取影响因子最高的期刊列表"""
    _load_journal_metrics()
    if not _journal_metrics:
        return []
    sorted_journals = sorted(
        _journal_metrics.values(),
        key=lambda x: x.get("jif", 0),
        reverse=True
    )
    return sorted_journals[:limit]


def ensure_dirs():
    ensure_workspace()
    DIGEST_DIR.mkdir(parents=True, exist_ok=True)
    for path in (PAPERS_DB, CONCEPTS_DB, EDGES_DB):
        path.touch(exist_ok=True)


def generate_paper_id(paper: dict) -> str:
    """生成论文唯一ID"""
    if paper.get("pmid"):
        return f"PMID:{paper['pmid']}"
    if paper.get("arxiv_id"):
        return f"ARXIV:{paper['arxiv_id']}"
    if paper.get("doi"):
        return f"DOI:{paper['doi']}"
    # biorxiv/medrxiv papers often have a DOI
    # fallback: title + first_author hash
    key = paper.get("title", "") + paper.get("first_author", "")
    return f"HASH:{hashlib.md5(key.encode()).hexdigest()[:12]}"


def paper_exists(paper_id: str) -> bool:
    """检查论文是否已索引"""
    if not PAPERS_DB.exists():
        return False
    with open(PAPERS_DB, "r", encoding="utf-8-sig") as f:
        for line in f:
            try:
                record = json.loads(line)
                if record.get("id") == paper_id:
                    return True
            except json.JSONDecodeError:
                continue
    return False


def add_paper(paper: dict) -> bool:
    """添加论文到知识图谱（去重）"""
    ensure_dirs()
    paper_id = generate_paper_id(paper)
    if paper_exists(paper_id):
        return False

    entry = {
        "id": paper_id,
        "title": paper.get("title", ""),
        "authors": paper.get("authors", []),
        "first_author": paper.get("first_author", ""),
        "last_author": paper.get("last_author", ""),
        "journal": paper.get("journal", ""),
        "year": paper.get("year", 0),
        "impact_factor": paper.get("impact_factor", 0.0),
        "impact_factor_5y": paper.get("impact_factor_5y", 0.0),
        "journal_quartile": paper.get("journal_quartile", ""),
        "journal_category": paper.get("journal_category", ""),
        "doi": paper.get("doi", ""),
        "pmid": paper.get("pmid", ""),
        "arxiv_id": paper.get("arxiv_id", ""),
        "abstract": paper.get("abstract", ""),
        "keywords": paper.get("keywords", []),
        "core_findings": paper.get("core_findings", ""),
        "core_question": paper.get("core_question", ""),
        "methods": paper.get("methods", []),
        # v3.0 深度分析字段
        "analysis_tier": paper.get("analysis_tier", "B"),  # S/A/B
        "claims": paper.get("claims", []),  # list of {claim, evidence, synthesis, strength, uncertain}
        "mechanism_cascade": paper.get("mechanism_cascade", ""),  # full text of mechanism cascade
        "hidden_axis": paper.get("hidden_axis", ""),  # hidden organizing axis
        "concept_innovation": paper.get("concept_innovation", ""),  # conceptual innovation
        "cross_refs": paper.get("cross_refs", []),  # list of related PMIDs/DOIs
        "evidence_strength": paper.get("evidence_strength", 0.0),  # average ★ score
        "entities": {
            "genes": paper.get("genes", []),
            "pathways": paper.get("pathways", []),
            "cell_types": paper.get("cell_types", []),
            "diseases": paper.get("diseases", []),
            "drugs": paper.get("drugs", []),
            "technologies": paper.get("technologies", []),
        },
        "relevance_score": paper.get("relevance_score", 0.0),
        "relevance_reason": paper.get("relevance_reason", ""),
        "source": paper.get("source", "unknown"),  # pubmed, arxiv, biorxiv, medrxiv, manual
        "fetched_at": datetime.now().isoformat(),
        "publication_date": paper.get("publication_date", ""),
    }

    with open(PAPERS_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    return True


def add_concept(concept: dict):
    """添加概念到概念库"""
    ensure_dirs()
    with open(CONCEPTS_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps({
            **concept,
            "added_at": datetime.now().isoformat()
        }, ensure_ascii=False) + "\n")


def add_edge(source_id: str, target_id: str, relation: str, metadata: dict = None):
    """添加关系边"""
    ensure_dirs()
    with open(EDGES_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps({
            "source": source_id,
            "target": target_id,
            "relation": relation,
            "metadata": metadata or {},
            "added_at": datetime.now().isoformat()
        }, ensure_ascii=False) + "\n")


def get_stats() -> dict:
    """获取知识图谱统计"""
    stats = {}
    for db_path, label in [(PAPERS_DB, "papers"), (CONCEPTS_DB, "concepts"), (EDGES_DB, "edges")]:
        if db_path.exists():
            with open(db_path, "r", encoding="utf-8-sig") as f:
                stats[label] = sum(1 for _ in f)
        else:
            stats[label] = 0

    # 按来源统计论文
    sources = {}
    if PAPERS_DB.exists():
        with open(PAPERS_DB, "r", encoding="utf-8-sig") as f:
            for line in f:
                try:
                    r = json.loads(line)
                    src = r.get("source", "unknown")
                    sources[src] = sources.get(src, 0) + 1
                except json.JSONDecodeError:
                    continue
    stats["by_source"] = sources

    # 最近更新
    stats["last_updated"] = datetime.now().isoformat()
    return stats


def write_daily_digest(content: str, date_str: str = None):
    """写入每日简报"""
    ensure_dirs()
    if date_str is None:
        date_str = date.today().isoformat()
    path = DIGEST_DIR / f"{date_str}.md"
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return str(path)


def get_recent_papers(days: int = 7) -> list:
    """获取最近 N 天的论文"""
    cutoff = datetime.now() - timedelta(days=days)
    papers = []
    if PAPERS_DB.exists():
        with open(PAPERS_DB, "r", encoding="utf-8-sig") as f:
            for line in f:
                try:
                    paper = json.loads(line)
                    fetched = datetime.fromisoformat(paper["fetched_at"])
                    if fetched >= cutoff:
                        papers.append(paper)
                except (json.JSONDecodeError, KeyError, ValueError):
                    continue
    return papers


def search_papers(query: str, field: str = "title") -> list:
    """简单搜索论文（关键词匹配）"""
    results = []
    if PAPERS_DB.exists():
        with open(PAPERS_DB, "r", encoding="utf-8-sig") as f:
            for line in f:
                try:
                    paper = json.loads(line)
                    text = ""
                    if field == "title":
                        text = paper.get("title", "")
                    elif field == "abstract":
                        text = paper.get("abstract", "")
                    elif field == "all":
                        text = json.dumps(paper, ensure_ascii=False)
                    
                    if query.lower() in text.lower():
                        results.append(paper)
                except json.JSONDecodeError:
                    continue
    return results


# ===== 用户研究领域定义 =====
# Customize this dict with your own research domains.
# The search-query-packs.json in assets/data/ provides pre-built query templates.
USER_RESEARCH_DOMAINS = {
    "primary": [],
    "secondary": [],
    "methods": [],
}


def get_pubmed_queries(date_str: str = None, domains: dict = None) -> list:
    """
    Generate PubMed search queries from research domains.

    Reads from USER_RESEARCH_DOMAINS or a passed domains dict.
    Format: {"primary": [...], "secondary": [...], "methods": [...]}
    Each entry generates one PubMed query filtered by date.
    """
    if domains is None:
        domains = USER_RESEARCH_DOMAINS
    
    if date_str is None:
        date_str = (date.today() - timedelta(days=1)).strftime("%Y/%m/%d")
    
    queries = []
    
    category_labels = {
        "primary": "核心",
        "secondary": "扩展",
        "methods": "方法",
    }
    
    for category, label in category_labels.items():
        for term in domains.get(category, []):
            queries.append({
                "label": f"{label}: {term}",
                "query": f'("{term}"[Title/Abstract]) AND ("{date_str}"[Date - Publication] : "{date_str}"[Date - Publication])',
            })
    
    return queries


# ===== 持久化辅助函数 =====

def load_existing_ids(db_path: Path) -> set:
    """加载已有论文ID集合，用于去重"""
    ids = set()
    if db_path.exists():
        with open(db_path, "r", encoding="utf-8-sig") as f:
            for line in f:
                try:
                    record = json.loads(line)
                    ids.add(record.get("id", ""))
                except json.JSONDecodeError:
                    continue
    return ids


def append_paper(paper: dict) -> bool:
    """追加论文到知识图谱（去重）"""
    ensure_dirs()
    paper_id = generate_paper_id(paper)
    existing = load_existing_ids(PAPERS_DB)
    if paper_id in existing:
        return False

    entry = dict(paper)
    entry["id"] = paper_id
    entry.setdefault("source", "unknown")
    entry.setdefault("fetched_at", datetime.now().isoformat())
    if "entities" not in entry:
        entry["entities"] = {
            "genes": paper.get("genes", []),
            "proteins": paper.get("proteins", []),
            "pathways": paper.get("pathways", []),
            "cell_types": paper.get("cell_types", []),
            "diseases": paper.get("diseases", []),
            "methods": paper.get("methods", []) or paper.get("technologies", []),
        }
    if entry.get("journal") and not entry.get("impact_factor"):
        entry = enrich_paper_if(entry)
    
    with open(PAPERS_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    return True


def append_concept(concept: dict) -> bool:
    """追加概念到概念库（去重）"""
    ensure_dirs()
    existing = load_existing_ids(CONCEPTS_DB)
    concept_id = concept.get("id", concept.get("name", ""))
    if concept_id in existing:
        return False
    
    with open(CONCEPTS_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps({
            "id": concept_id,
            **concept,
            "added_at": datetime.now().isoformat()
        }, ensure_ascii=False) + "\n")
    return True


def append_edge(source_id: str, target_id: str, relation: str, metadata: dict = None):
    """追加关系边（去重）"""
    ensure_dirs()
    edges = []
    if EDGES_DB.exists():
        with open(EDGES_DB, "r", encoding="utf-8-sig") as f:
            for line in f:
                try:
                    e = json.loads(line)
                    if e.get("source") == source_id and e.get("target") == target_id:
                        return False  # already exists
                except json.JSONDecodeError:
                    continue
    
    with open(EDGES_DB, "a", encoding="utf-8") as f:
        f.write(json.dumps({
            "source": source_id,
            "target": target_id,
            "relation": relation,
            "metadata": metadata or {},
            "added_at": datetime.now().isoformat()
        }, ensure_ascii=False) + "\n")
    return True


# 兼容旧接口
add_paper = append_paper
add_concept = append_concept
add_edge = append_edge


# ===== CLI =====
if __name__ == "__main__":
    import sys
    cmd = sys.argv[1] if len(sys.argv) > 1 else "stats"
    
    if cmd == "stats":
        stats = get_stats()
        print(json.dumps(stats, indent=2, ensure_ascii=False))
    elif cmd == "recent":
        days = int(sys.argv[2]) if len(sys.argv) > 2 else 7
        papers = get_recent_papers(days)
        print(f"最近 {days} 天: {len(papers)} 篇论文")
        for p in papers:
            print(f"  [{p['id']}] {p['title'][:80]}")
    elif cmd == "search":
        query = sys.argv[2] if len(sys.argv) > 2 else ""
        papers = search_papers(query)
        print(f"搜索 '{query}': {len(papers)} 篇")
        for p in papers[:20]:
            print(f"  [{p['id']}] {p['title'][:80]}")
    elif cmd == "queries":
        queries = get_pubmed_queries()
        for q in queries:
            print(f"\n# {q['label']}")
            print(q['query'])
    elif cmd == "network":
        import subprocess
        # Run build_network.py from the scripts directory
        scripts_dir = Path(__file__).resolve().parent
        subprocess.run([sys.executable, str(scripts_dir / "build_network.py")])
    else:
        print(f"Unknown command: {cmd}")
        print("Usage: python kg_core.py [stats|recent|search|queries|network]")
