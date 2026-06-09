"""
知识图谱关联边生成器 v3.1 — Bioconductor 驱动的语义关联
====================================================
五层生物学关联逻辑（全部倒排索引，无 O(N²) 遍历）：

  策略 1  显式引用    tier7_cross_refs + related_pmids（排除 same_journal/issue）
  策略 2  共享分子     org.Hs+Mm.eg.db 73k 基因符号 + KEGG 通路（两段式匹配）
  策略 2.5 文本重叠    bag-of-words 倒排索引，≥4 共享关键词
  策略 3  同病+同法    疾病标签 × 方法标签交叉
  策略 4  隐藏轴      Tier5 深层范式共鸣
  策略 5  概念节点    concepts.db → 论文

Data dependencies live under $LITERATURE_KG_ROOT/data:
  bioc_genes.json and kegg_pathways.json
生成方式: wsl -d Ubuntu-20.04 conda activate sc_spatial_env
          Rscript -e 'library(org.Hs.eg.db); library(org.Mm.eg.db); ...'
"""

import json, re
from pathlib import Path
from collections import defaultdict, Counter
from itertools import combinations

from workspace_paths import DATA_ASSETS_DIR, KG_ROOT

DATA_DIR = KG_ROOT / "data"
PAPERS_DB = KG_ROOT / "papers.db"
CONCEPTS_DB = KG_ROOT / "concepts.db"
EDGES_DB = KG_ROOT / "edges.db"

# ═══════════════════════════════════════════════
#  Bioconductor 词表加载（懒加载，单例）
# ═══════════════════════════════════════════════

def _load_bioc_genes():
    """org.Hs.eg.db + org.Mm.eg.db 导出的人+鼠基因符号集合"""
    path = DATA_DIR / "bioc_genes.json"
    if not path.exists():
        path = DATA_ASSETS_DIR / "bioc_genes.json"
    if path.exists():
        with open(path, 'r', encoding='utf-8-sig') as f:
            return set(g.upper() for g in json.load(f))
    return set()

def _load_kegg_pathways():
    """KEGGREST + GO.db 导出的通路/GO术语（KEGG 297 + GO BP 25,647）"""
    path = DATA_DIR / "kegg_pathways.json"
    if not path.exists():
        path = DATA_ASSETS_DIR / "kegg_pathways.json"
    if path.exists():
        with open(path, 'r', encoding='utf-8-sig') as f:
            return json.load(f)
    return []

_BIOC_GENES = None
_KEGG_PATHWAYS = None

def bioc_genes():
    global _BIOC_GENES
    if _BIOC_GENES is None:
        _BIOC_GENES = _load_bioc_genes()
    return _BIOC_GENES

def kegg_pathways():
    global _KEGG_PATHWAYS
    if _KEGG_PATHWAYS is None:
        _KEGG_PATHWAYS = _load_kegg_pathways()
    return _KEGG_PATHWAYS

# ═══════════════════════════════════════════════
#  高频分子补充（Bioconductor ASCII 别名覆盖不到的希腊字母/特殊符号变体）
# ═══════════════════════════════════════════════

HIGH_CONFIDENCE_EXTRAS = {
    'NF-κB', 'TGF-β', 'TNF-α', 'TNF', 'IL-1β', 'IL-6', 'IL-10', 'IL-17',
    'IFN-α', 'IFN-β', 'IFN-γ', 'Tau', 'Aβ', 'α-synuclein',
    'TREM2', 'ApoE', 'C1q', 'C3',
    'CD40', 'CD40L', 'CTLA-4', 'PD-1', 'PD-L1', 'LAG3', 'TIM3', 'TIGIT',
    'CD19', 'CD20', 'BCMA', 'GPRC5D', 'CD3', 'CD4', 'CD8',
    'NFAT5', 'PPARγ', 'p53', 'Axl', 'Gas6', 'MerTK',
    'Piezo1', 'Piezo2', 'SIRPα', 'GLP-1', 'GLP-1R', 'GIP', 'GIPR',
    'MG53', 'GPX4', 'FSP1', 'DHODH', 'SLC7A11',
    'p62', 'SQSTM1', 'BECN1', 'ATG5', 'ATG7', 'ATG12', 'ULK1',
}

# ═══════════════════════════════════════════════
#  正则
# ═══════════════════════════════════════════════

RE_GENE_ABBREV = re.compile(r'\b([A-Z]{2,6}\d{0,3}[A-Z]?)\b')
RE_TITLECASE_GENE = re.compile(r'\b([A-Z][a-z]{2,7}\d{0,3}[a-z]?)\b')  # 小鼠: Tnf, Il6, Trp53
RE_MIRNA = re.compile(r'\b(miR-\d+[a-z]?-\d+p?)\b')
RE_PATHWAY = re.compile(
    r'(Wnt|Hedgehog|Notch|TGF-?[βB]|NF-?[κk]B|MAPK|PI3K[/-]Akt|JAK[/-]STAT|mTOR|AMPK|Hippo|'
    r'RAS|RAF|RTK|GPCR|JNK|p38|ERK|'
    r'cGAS-STING|inflammasome|complement|coagulation|fibrinolysis|'
    r'ferroptosis|autophagy|apoptosis|pyroptosis|necroptosis|'
    r'glycolysis|OXPHOS|PPP|TCA)\s+(pathway|signaling|axis|signalling|cascade)',
    re.IGNORECASE)

STOP_MOLECULES = {
    'PHASE','TRIAL','THERAPY','TREATMENT','STUDY','PATIENT','GROUP','RESULT','DATA',
    'USING','BASED','ROLE','EFFECT','LEVEL','YEAR','WEEK','MONTH','DAY',
    'FOUND','SHOW','ALSO','WELL','MAY','ONE','TWO','NEW','KEY','MAIN','HIGH','LOW',
    'CELL','HUMAN','MOUSE','MODEL','TYPE','FORM','LIKE','VIA','DUE','AGE','SEX',
    'GENE','PROTEIN','TARGET','DRUG','DOSE','SITE','LOSS','GAIN','CHANGE',
    'INCREASE','DECREASE','METHOD','RESULT','CONCLUSION','OBJECTIVE','BACKGROUND',
    'THESE','THOSE','THEIR','THEY','WERE','FROM','HAVE','BEEN','MORE','LESS',
    'ANALYSIS','ANALYSES','CLINICAL','MEDICAL','MAJOR','MINOR','FIRST','LAST',
    'RECENT','PREVIOUS','EARLY','LATE','PRIMARY','SECONDARY','SINGLE','DOUBLE',
    'PLACEBO','CONTROL','RCT','MRI','CT','PET','RNA','DNA','MRNA','LNCRNA',
    'ASSOCIATED','INCLUDING','POTENTIAL','FUNCTIONAL','DIFFERENT','SIGNIFICANT',
    'IMPORTANT','CRITICAL','ESSENTIAL','SPECIFIC','ADDITIONAL',
    'CURRENTLY','CURRENT','RECENTLY','FREQUENTLY','TYPICALLY','GENERALLY',
    'NATURE','SCIENCE','LANCET','JAMA','BMJ','NEJM','TRUE','FALSE','YES','NO',
    # 非基因缩写
    'USA','UK','FDA','EMA','NIH','WHO','CDC','AIM','AIMS','ART','BMI','BP','CI',
    'COX','CRP','CT','DOI','ECG','ELISA','FACS','FISH','GCP','HIV','HLA','IHC',
    'ISH','IV','MAP','MRI','NGS','PCR','PK','PD','QOL','RCT','ROC','RT','SEM',
    'SNP','TMA','TMB','WES','WGS','WT','KO','KD','OE','IP','WB','IF',
    'COVID','SARS','MERS','HIV','HBV','HCV','HPV','EBV','CMV',
    'AUC','NPV','PPV','MCC','MSE','MAE','RMSE','OR','HR','RR','SD','SE',
}

DISEASE_PATTERNS = [
    r'\b(Parkinson|Alzheimer|cancer|tumor|diabetes|obesity|stroke|CKD|COPD|SLE|ALS|IPF|HCC|GBM|TB|MS|RA|IBD|CVD|CAD|HF|ESRD|NAFLD|NASH|PCOS|ADHD|PTSD|OCD|MDD)\b',
    r'\b(fibrosis|inflammation|neurodegeneration|metastasis|resistance|autoimmune)\b',
]

METHOD_PATTERNS = [
    r'\b(CRISPR|RNA-seq|scRNA-seq|ATAC-seq|ChIP-seq|GWAS|MR|cryo-EM|organoid|iPSC|PDX)\b',
    r'\b(single.cell|spatial.transcriptom|proteomics|metabolomics|machine.learning|deep.learning)\b',
]

# ═══════════════════════════════════════════════
#  实体提取
# ═══════════════════════════════════════════════

def extract_entities(text):
    """两段式基因匹配：正则提取候选 → 查 Bioconductor 集合（O(候选数)）"""
    entities = {'molecules': set(), 'diseases': set(), 'methods': set()}
    genes = bioc_genes()

    # 高频额外分子（希腊字母/特殊符号变体）
    utext = text.upper()
    for mol in HIGH_CONFIDENCE_EXTRAS:
        if mol.upper() in utext:
            entities['molecules'].add(mol.upper())

    # 全大写基因 (TNF, EGFR, BCL2)
    for m in RE_GENE_ABBREV.findall(text):
        m_up = m.upper()
        if m_up not in STOP_MOLECULES and m_up in genes:
            entities['molecules'].add(m_up)

    # TitleCase 小鼠基因 (Tnf, Il6, Smad2)
    for m in RE_TITLECASE_GENE.findall(text):
        m_up = m.upper()
        if m_up not in STOP_MOLECULES and m_up in genes:
            entities['molecules'].add(m_up)

    # miRNA
    for m in RE_MIRNA.findall(text):
        entities['molecules'].add(m.upper())

    # 通路名
    for m in RE_PATHWAY.findall(text):
        full = ' '.join(m) if isinstance(m, tuple) else m
        clean = full.strip().upper().replace(' PATHWAY','').replace(' SIGNALING','').replace(' AXIS','').replace(' CASCADE','')
        if len(clean) >= 4:
            entities['molecules'].add(clean)

    # KEGG + GO 通路匹配（子串匹配，预 lowered 文本）
    ltext = text.lower()
    for kp in kegg_pathways():
        if kp.lower() in ltext:
            entities['molecules'].add(kp.upper())

    # 疾病 / 方法
    for pat in DISEASE_PATTERNS:
        for m in re.findall(pat, text, re.IGNORECASE):
            entities['diseases'].add(m.lower())
    for pat in METHOD_PATTERNS:
        for m in re.findall(pat, text, re.IGNORECASE):
            entities['methods'].add(m.lower())

    return entities

# ═══════════════════════════════════════════════
#  工具函数
# ═══════════════════════════════════════════════

def analysis_block(paper):
    value = paper.get("analysis")
    return value if isinstance(value, dict) else {}


def entity_block(paper):
    value = paper.get("entities")
    return value if isinstance(value, dict) else paper


def core_question(paper):
    analysis = analysis_block(paper)
    return (
        paper.get("tier2_core_question")
        or paper.get("core_question")
        or analysis.get("core_question")
        or ""
    )


def subquestions(paper):
    analysis = analysis_block(paper)
    return paper.get("tier2_subquestions") or analysis.get("subquestions") or []


def mechanism_block(paper):
    analysis = analysis_block(paper)
    return (
        paper.get("tier4_mechanism_cascade")
        or paper.get("mechanism_cascade")
        or analysis.get("mechanism")
        or {}
    )


def mechanism_steps(paper):
    mechanism = mechanism_block(paper)
    if isinstance(mechanism, dict):
        return mechanism.get("cascade") or mechanism.get("steps") or []
    return [mechanism] if mechanism else []


def hidden_axes(paper):
    analysis = analysis_block(paper)
    return (
        paper.get("tier5_hidden_axis")
        or paper.get("hidden_axis")
        or analysis.get("hidden_axes")
        or []
    )


def cross_references(paper):
    analysis = analysis_block(paper)
    return (
        paper.get("tier7_cross_refs")
        or paper.get("cross_refs")
        or analysis.get("cross_references")
        or []
    )


def curated_entities(paper):
    values = {"molecules": set(), "diseases": set(), "methods": set()}
    container = entity_block(paper)
    for field in ("genes", "proteins", "pathways"):
        values["molecules"].update(
            str(value).upper() for value in (container.get(field) or []) if value
        )
    values["diseases"].update(
        str(value).lower() for value in (container.get("diseases") or []) if value
    )
    for field in ("methods", "technologies"):
        values["methods"].update(
            str(value).lower() for value in (container.get(field) or []) if value
        )
    return values


def paper_text(paper):
    """聚合论文文本——优先读取全文缓存"""
    pid = paper.get('id', '')
    
    # 缓存查找
    safe = re.sub(r'[^a-zA-Z0-9_.-]', '_', pid)
    cache_file = KG_ROOT / "fulltext_cache" / f"{safe}.txt"
    if cache_file.exists():
        text = cache_file.read_text(encoding='utf-8-sig', errors='replace')
        text = re.sub(r'^#.*?\n\n', '', text, flags=re.DOTALL)
        if len(text) > 200:
            return text
    
    # 回退：数据库字段
    parts = []
    for k in ['title', 'abstract', 'core_findings']:
        v = paper.get(k, '')
        if v and not str(v).startswith('非研究'):
            parts.append(str(v))
    question = core_question(paper)
    if question:
        parts.append(str(question))
    analysis = analysis_block(paper)
    for chain in analysis.get("ces_chains", []) or []:
        if isinstance(chain, dict):
            parts.extend(
                str(chain.get(key, ""))
                for key in ("claim", "evidence", "synthesis")
                if chain.get(key)
            )
    parts.extend(str(step) for step in mechanism_steps(paper) if step)
    ha = hidden_axes(paper)
    if isinstance(ha, list):
        for item in ha:
            if isinstance(item, dict):
                parts.append(str(item.get('observation', '')))
                parts.append(str(item.get('interpretation', '')))
    elif isinstance(ha, dict):
        for k in ('observation_1', 'observation_2', 'observation_3', 'interpretation_1'):
            if ha.get(k):
                parts.append(str(ha[k]))
    return ' '.join(parts)

def load_papers():
    papers = []
    if not PAPERS_DB.exists():
        return papers
    with open(PAPERS_DB, 'r', encoding='utf-8-sig') as f:
        for line in f:
            if not line.strip(): continue
            try: papers.append(json.loads(line))
            except: continue
    return papers

def pick_title(p):
    return (p.get('title') or '')[:60].strip()

def pick_topic(p, max_len=80):
    cq = core_question(p)
    if cq:
        return re.sub(r'[？?].*', '', cq)[:max_len].strip()
    cf = p.get('core_findings', '')
    if cf and not cf.startswith('非研究'):
        return re.sub(r'[。.]', '', cf)[:max_len].strip()
    return pick_title(p)

# ═══════════════════════════════════════════════
#  主逻辑
# ═══════════════════════════════════════════════

def main():
    print("加载 Bioconductor 词表 ...")
    bg = bioc_genes()
    kp = kegg_pathways()
    print(f"  基因符号: {len(bg):,}, 通路(KEGG+GO): {len(kp):,}")

    print("加载论文 ...")
    papers = load_papers()
    print(f"  已加载 {len(papers)} 篇")
    papers = [p for p in papers if p.get("id")]
    id_map = {p['id']: p for p in papers}

    print("计算实体指纹 ...")
    fingerprints = {}
    for p in papers:
        fingerprint = extract_entities(paper_text(p))
        curated = curated_entities(p)
        for field in fingerprint:
            fingerprint[field].update(curated[field])
        fingerprints[p['id']] = fingerprint
    print(f"  完成 {len(fingerprints)} 篇")

    # 加载已有边
    existing_pairs = set()
    existing_edges = []
    if EDGES_DB.exists():
        with open(EDGES_DB, 'r', encoding='utf-8-sig') as f:
            for line in f:
                if not line.strip(): continue
                try:
                    e = json.loads(line)
                    existing_pairs.add((e['source'], e['target']))
                    existing_edges.append(e)
                except: pass

    new_edges = []

    def add_edge(src, tgt, rel, desc):
        key = (src, tgt)
        if key in existing_pairs:
            return
        existing_pairs.add(key)
        new_edges.append({'source': src, 'target': tgt, 'relation': rel, 'description': desc})

    # ── 策略 1: 显式引用 ──
    print("策略 1: 显式引用 ...")
    NON_BIO_RELS = {'same_journal', 'same_issue', 'same_author'}
    for p in papers:
        pid = p['id']
        for ref in cross_references(p):
            if not isinstance(ref, dict):
                continue
            rel = ref.get('relation', 'cites')
            if rel in NON_BIO_RELS:
                continue
            tgt_ref = ref.get('ref_id', ref.get('pmid', ''))
            if tgt_ref:
                tgt_id = tgt_ref if tgt_ref.startswith('PMID:') or tgt_ref.startswith('BIORXIV:') or tgt_ref.startswith('ARXIV:') else f'PMID:{tgt_ref}'
                if tgt_id in id_map:
                    desc = ref.get('description', ref.get('desc', '')) or f'{pick_title(p)} <-> {pick_title(id_map[tgt_id])}'
                    add_edge(pid, tgt_id, rel, desc)
        for rp in p.get('related_pmids', []):
            tgt_id = f'PMID:{rp}' if not str(rp).startswith('PMID:') else str(rp)
            if tgt_id in id_map and tgt_id != pid:
                add_edge(pid, tgt_id, 'cites',
                         f'「{pick_topic(p, 60)}」引用「{pick_topic(id_map[tgt_id], 60)}」')

    # ── 策略 2: 共享分子实体（倒排索引） ──
    print("策略 2: 共享分子 ...")
    deep_papers = [p for p in papers if mechanism_steps(p)]

    mol_index = defaultdict(set)
    for p in deep_papers:
        for m in fingerprints[p['id']]['molecules']:
            if m not in STOP_MOLECULES and len(m) >= 3:
                mol_index[m].add(p['id'])

    MAX_PER_MOL = 15
    pair_shared = defaultdict(lambda: {'mols': set()})
    for mol, pids in mol_index.items():
        if len(pids) < 2: continue
        pid_list = list(pids)[:MAX_PER_MOL]
        for i in range(len(pid_list)):
            for j in range(i+1, len(pid_list)):
                key = tuple(sorted([pid_list[i], pid_list[j]]))
                pair_shared[key]['mols'].add(mol)

    for (pid1, pid2), data in pair_shared.items():
        shared = data['mols']
        if len(shared) >= 2:
            shared_str = ', '.join(sorted(shared)[:5])
            add_edge(pid1, pid2, 'shares_molecules',
                     f'共享分子「{shared_str}」: 「{pick_topic(id_map[pid1], 50)}」|「{pick_topic(id_map[pid2], 50)}」')

    print(f"  deep_papers={len(deep_papers)}, mol_index={len(mol_index)}, pair_shared={len(pair_shared)}")

    # ── 策略 2.5: 文本主题重叠（倒排索引） ──
    print("策略 2.5: 文本重叠 ...")
    STOP_WORDS = {
        'the','a','an','is','are','was','were','be','been','being','have','has','had',
        'do','does','did','will','would','could','should','may','might','can','shall',
        'to','of','in','for','on','with','at','by','from','as','into','through','during',
        'and','or','not','but','if','while','than','nor','so','yet','both','either',
        'this','that','these','those','it','its','they','them','their','we','our','us',
        'which','who','whom','what','when','where','how','all','each','every','some','any',
        'no','none','other','another','such','only','own','same','more','most','also',
        'very','just','now','then','here','there','up','out','off','new','first',
        'using','based','found','show','shown','reported','identified','demonstrate',
        'associated','including','increased','decreased','significantly','compared',
        'however','therefore','thus','furthermore','recent','previous','current',
        'study','studies','research','findings','results','data','analysis','method',
        'approach','model','role','effect','level','group','patient','cell','human',
        'mouse','treatment','clinical','therapeutic','potential','important','key',
    }

    keywords_by_paper = {}
    for p in papers:
        text = f"{p.get('title','')} {p.get('core_findings','')} {core_question(p)}".lower()
        kw = set(re.findall(r'[a-z]{4,}', text)) - STOP_WORDS
        if kw:
            keywords_by_paper[p['id']] = (kw, p)

    word_index = defaultdict(set)
    for pid, (kws, _) in keywords_by_paper.items():
        for w in kws:
            word_index[w].add(pid)

    text_pairs = defaultdict(lambda: {'words': set()})
    for word, pids in word_index.items():
        if len(pids) < 2 or len(pids) > 100: continue
        pid_list = list(pids)[:30]
        for i in range(len(pid_list)):
            for j in range(i+1, len(pid_list)):
                key = tuple(sorted([pid_list[i], pid_list[j]]))
                text_pairs[key]['words'].add(word)

    for (pid1, pid2), data in text_pairs.items():
        shared = data['words']
        if len(shared) >= 4:
            shared_str = ', '.join(sorted(shared)[:6])
            add_edge(pid1, pid2, 'shares_topic',
                     f'主题相似「{shared_str}」: 「{pick_topic(id_map[pid1], 40)}」⇄「{pick_topic(id_map[pid2], 40)}」')

    print(f"  text_pairs={len(text_pairs)}, word_index={len(word_index)}")

    # ── 策略 3: 同病+同法 ──
    print("策略 3: 同病+同法 ...")
    disease_index = defaultdict(set)
    method_index = defaultdict(set)
    for p in deep_papers:
        fp = fingerprints[p['id']]
        for d in fp['diseases']:
            disease_index[d].add(p['id'])
        for m in fp['methods']:
            method_index[m].add(p['id'])

    disease_pairs = defaultdict(int)
    for disease, pids in disease_index.items():
        if len(pids) < 2: continue
        pid_list = list(pids)[:20]
        for i in range(len(pid_list)):
            for j in range(i+1, len(pid_list)):
                disease_pairs[tuple(sorted([pid_list[i], pid_list[j]]))] += 1

    for meth, pids in method_index.items():
        if len(pids) < 2: continue
        pid_list = list(pids)[:20]
        for i in range(len(pid_list)):
            for j in range(i+1, len(pid_list)):
                key = tuple(sorted([pid_list[i], pid_list[j]]))
                if key in disease_pairs:
                    sd = fingerprints[pid_list[i]]['diseases'] & fingerprints[pid_list[j]]['diseases']
                    sm = fingerprints[pid_list[i]]['methods'] & fingerprints[pid_list[j]]['methods']
                    add_edge(key[0], key[1], 'shares_disease_method',
                             f'同病「{", ".join(sorted(sd)[:3])}」+同法「{", ".join(sorted(sm)[:3])}」: '
                             f'「{pick_topic(id_map[key[0]], 40)}」⇄「{pick_topic(id_map[key[1]], 40)}」')

    # ── 策略 4: Tier5 隐藏轴共鸣 ──
    print("策略 4: 隐藏轴 ...")
    MEANINGFUL = {'paradigm','principle','bias','model','framework','hidden','axis',
                  'survivor','identity','recognition'}
    for p1, p2 in combinations(deep_papers[:200], 2):
        ha1, ha2 = hidden_axes(p1), hidden_axes(p2)
        if not ha1 or not ha2: continue
        def _ha_text(ha):
            if isinstance(ha, list):
                return ' '.join(str(v) for item in ha if isinstance(item, dict) for v in item.values() if isinstance(v, str))
            if isinstance(ha, dict):
                return ' '.join(str(v) for v in ha.values() if isinstance(v, str))
            return str(ha)
        t1 = _ha_text(ha1)
        t2 = _ha_text(ha2)
        if set(t1.lower().split()) & set(t2.lower().split()) & MEANINGFUL:
            add_edge(p1['id'], p2['id'], 'shares_paradigm',
                     f'隐藏轴共鸣: 「{pick_title(p1)}」⇄「{pick_title(p2)}」')

    # ── 策略 5: 概念节点 ──
    print("策略 5: 概念节点 ...")
    if CONCEPTS_DB.exists():
        with open(CONCEPTS_DB, 'r', encoding='utf-8-sig') as f:
            for line in f:
                if not line.strip(): continue
                try:
                    c = json.loads(line)
                    source_papers = c.get('source_papers') or []
                    if c.get('source_paper'):
                        source_papers = [c['source_paper'], *source_papers]
                    for source_paper in dict.fromkeys(source_papers):
                        if source_paper and c.get('id'):
                            add_edge(source_paper, c['id'], 'defines_concept',
                                     f'定义概念「{c.get("name", c["id"])}」: {c.get("definition", "")[:120]}')
                except: pass

    # ── 写入 ──
    print("写入 edges.db ...")
    with open(EDGES_DB, 'w', encoding='utf-8') as f:
        for e in existing_edges + new_edges:
            f.write(json.dumps(e, ensure_ascii=False) + '\n')

    # ── 统计报告 ──
    total = len(existing_edges) + len(new_edges)
    rels = Counter(e['relation'] for e in new_edges)
    print(f"\n已有边: {len(existing_edges)}  +  新增: {len(new_edges)}  =  总计: {total}")
    print("新增边类型:")
    for r, c in rels.most_common():
        print(f"  {r}: {c}")
    if new_edges:
        print("\n新知边样例:")
        for e in new_edges[:6]:
            print(f"  [{e['relation']}] {e['description'][:150]}")


if __name__ == "__main__":
    main()
