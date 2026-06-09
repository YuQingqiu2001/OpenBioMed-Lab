from pathlib import Path
from collections import Counter, defaultdict
from datetime import date
import json
import re
import time
import hashlib
import urllib.request

from workspace_paths import KG_ROOT

ROOT = KG_ROOT
REPORT = {
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'root': str(ROOT),
    'inventory': {},
    'forbidden_residue': {},
    'ports': {},
    'db': {},
    's_tier': {},
    'edges': {},
    'concepts': {},
    'fulltext_cache': {},
    'issues': [],
    'warnings': [],
}


def fmt(n):
    n = float(n)
    for u in ['B', 'KB', 'MB', 'GB', 'TB']:
        if n < 1024 or u == 'TB':
            return f'{n:.1f}{u}' if u != 'B' else f'{int(n)}B'
        n /= 1024


def add_issue(sev, msg, data=None):
    item = {'severity': sev, 'message': msg}
    if data is not None:
        item['data'] = data
    if sev in ['critical', 'error']:
        REPORT['issues'].append(item)
    else:
        REPORT['warnings'].append(item)


def analysis_view(paper):
    analysis = paper.get('analysis')
    analysis = analysis if isinstance(analysis, dict) else {}
    return {
        'subquestions': paper.get('tier2_subquestions') or analysis.get('subquestions') or [],
        'ces_chains': paper.get('tier3_ces_chains') or analysis.get('ces_chains') or [],
        'mechanism': (
            paper.get('tier4_mechanism_cascade')
            or paper.get('mechanism_cascade')
            or analysis.get('mechanism')
            or {}
        ),
        'hidden_axes': (
            paper.get('tier5_hidden_axis')
            or paper.get('hidden_axis')
            or analysis.get('hidden_axes')
            or []
        ),
        'conceptual_contribution': (
            paper.get('tier6_concept_innovation')
            or paper.get('concept_innovation')
            or analysis.get('conceptual_contribution')
            or {}
        ),
        'cross_references': (
            paper.get('tier7_cross_refs')
            or paper.get('cross_refs')
            or analysis.get('cross_references')
            or []
        ),
    }


def has_content(value):
    if isinstance(value, dict):
        return any(has_content(item) for item in value.values())
    if isinstance(value, (list, tuple, set)):
        return any(has_content(item) for item in value)
    return bool(str(value).strip()) if value is not None else False


def analysis_tier(paper):
    explicit = paper.get('analysis_tier')
    if explicit:
        return explicit
    return 'S' if any(has_content(value) for value in analysis_view(paper).values()) else '<missing>'


def scan_inventory():
    files = []
    for p in ROOT.rglob('*'):
        try:
            if p.is_file():
                st = p.stat()
                files.append({
                    'path': p,
                    'rel': p.relative_to(ROOT).as_posix(),
                    'size': st.st_size,
                    'mtime': st.st_mtime,
                    'suffix': p.suffix.lower(),
                })
        except Exception as e:
            add_issue('warning', f'stat failed: {p}', repr(e))
    REPORT['inventory']['total_files'] = len(files)
    total = sum(f['size'] for f in files)
    REPORT['inventory']['total_size_bytes'] = total
    REPORT['inventory']['total_size_human'] = fmt(total)
    sizes = defaultdict(int)
    counts = defaultdict(int)
    for f in files:
        top = f['rel'].split('/')[0]
        sizes[top] += f['size']
        counts[top] += 1
    REPORT['inventory']['top_level'] = [
        {'name': k, 'files': counts[k], 'size': fmt(sizes[k])}
        for k in sorted(sizes, key=sizes.get, reverse=True)
    ]
    return files


def check_forbidden_and_ports():
    forbidden = [
        'chrome_cdp_profile', '.chrome_biorxiv_profile', '.chrome_cdp_profile',
        '.biorxiv_cookie', 'biorxiv_api_test', 'biorxiv_api_test2',
        'biorxiv_api_test_doi', 'biorxiv_browser_test',
    ]
    for rel in forbidden:
        exists = (ROOT / rel).exists()
        REPORT['forbidden_residue'][rel] = exists
        if exists:
            add_issue('warning', f'Forbidden residue still exists: {rel}')
    for port in [9222, 9223]:
        try:
            with urllib.request.urlopen(f'http://127.0.0.1:{port}/json/version', timeout=1) as r:
                data = json.loads(r.read().decode('utf-8', 'replace'))
            REPORT['ports'][str(port)] = {
                'open': True,
                'browser': data.get('Browser'),
                'ua': data.get('User-Agent'),
            }
            add_issue('warning', f'CDP port {port} is still open', REPORT['ports'][str(port)])
        except Exception:
            REPORT['ports'][str(port)] = {'open': False}


def load_db(name):
    p = ROOT / name
    rows = []
    bad = []
    if not p.exists():
        add_issue('critical', f'Missing DB: {name}')
    else:
        with p.open('r', encoding='utf-8-sig') as f:
            for i, line in enumerate(f, 1):
                if not line.strip():
                    continue
                try:
                    rows.append(json.loads(line))
                except Exception as e:
                    bad.append({'line': i, 'error': repr(e), 'preview': line[:200]})
    REPORT['db'][name] = {
        'rows': len(rows),
        'bad_lines': len(bad),
        'bad_examples': bad[:5],
        'size': fmt(p.stat().st_size) if p.exists() else None,
    }
    if bad:
        add_issue('critical', f'{name} has bad NDJSON lines', bad[:3])
    return rows


def audit_papers(papers):
    ids = [p.get('id') for p in papers]
    id_counter = Counter(ids)
    dup_ids = [{'id': k, 'count': v} for k, v in id_counter.items() if k and v > 1]
    missing_ids = sum(1 for x in ids if not x)
    required = ['id', 'title', 'source']
    missing_req = []
    for i, p in enumerate(papers):
        miss = [k for k in required if not p.get(k)]
        if miss:
            missing_req.append({'idx': i, 'id': p.get('id'), 'missing': miss})
    REPORT['db']['papers_integrity'] = {
        'unique_ids': len([k for k in id_counter if k]),
        'missing_ids': missing_ids,
        'duplicate_ids': dup_ids[:20],
        'missing_required_examples': missing_req[:20],
        'source_counts': dict(Counter(p.get('source', '<missing>') for p in papers)),
        'analysis_tier_counts': dict(Counter(analysis_tier(p) for p in papers)),
    }
    if missing_ids or dup_ids:
        add_issue('error', 'Paper ID problems', {'missing_ids': missing_ids, 'duplicate_ids': dup_ids[:10]})
    if missing_req:
        add_issue('warning', 'Some papers missing basic required fields', missing_req[:10])


def audit_s_tier(papers):
    s_papers = [p for p in papers if analysis_tier(p) == 'S']
    empty_s = []
    weak_s = []
    counts = Counter()
    for p in s_papers:
        pid = p.get('id')
        view = analysis_view(p)
        t2 = view['subquestions']
        t3 = view['ces_chains']
        t4 = view['mechanism']
        t5 = view['hidden_axes']
        t6 = view['conceptual_contribution']
        t7 = view['cross_references']
        evidence_ok = 0
        if isinstance(t3, list):
            for c in t3:
                ev = str(c.get('evidence', '') if isinstance(c, dict) else '')
                if len(ev.strip()) > 20:
                    evidence_ok += 1
        rec = {
            'id': pid,
            'title': p.get('title', '')[:120],
            'T2': len(t2) if isinstance(t2, list) else -1,
            'T3': len(t3) if isinstance(t3, list) else -1,
            'T3_evidence_ok': evidence_ok,
            'T4': has_content(t4),
            'T5': has_content(t5),
            'T6': has_content(t6),
            'T7': len(t7) if isinstance(t7, list) else -1,
        }
        if rec['T2'] == 0 and rec['T3'] == 0 and not rec['T4']:
            empty_s.append(rec)
        elif not (rec['T2'] >= 5 and rec['T3'] >= 5 and rec['T3_evidence_ok'] >= 5 and rec['T4'] and rec['T5'] and rec['T6'] and rec['T7'] >= 5):
            weak_s.append(rec)
        counts['T2>=5'] += rec['T2'] >= 5
        counts['T3>=5'] += rec['T3'] >= 5
        counts['T3_ev>=5'] += rec['T3_evidence_ok'] >= 5
        counts['T4'] += rec['T4']
        counts['T5'] += rec['T5']
        counts['T6'] += rec['T6']
        counts['T7>=5'] += rec['T7'] >= 5
    REPORT['s_tier'] = {
        'total_S': len(s_papers),
        'field_pass_counts': dict(counts),
        'empty_s_count': len(empty_s),
        'weak_s_count': len(weak_s),
        'empty_examples': empty_s[:20],
        'weak_examples': weak_s[:30],
    }
    if empty_s:
        add_issue('error', 'Empty-shell S-tier papers detected', empty_s[:10])
    if weak_s:
        add_issue('warning', 'S-tier papers not meeting strict v4.0 completeness', {'count': len(weak_s), 'examples': weak_s[:10]})


def audit_concepts(concepts):
    concept_ids = [c.get('id') for c in concepts]
    concept_dups = [{'id': k, 'count': v} for k, v in Counter(concept_ids).items() if k and v > 1]
    concept_missing = [c for c in concepts if not c.get('id') or not c.get('name')]
    REPORT['concepts'] = {
        'rows': len(concepts),
        'duplicate_ids': concept_dups,
        'missing_id_or_name': concept_missing[:10],
    }
    if concept_dups or concept_missing:
        add_issue('warning', 'Concept problems', REPORT['concepts'])


def audit_edges(edges, papers, concepts):
    paper_ids = {p.get('id') for p in papers if p.get('id')}
    concept_ids = {c.get('id') for c in concepts if c.get('id')}
    node_ids = paper_ids | concept_ids
    illegal_relations = {'same_journal', 'same_issue', 'same_author', 'shares_keyword'}
    edge_keys = Counter()
    illegal = []
    no_desc = []
    orphan = []
    selfloops = []
    relations = Counter()
    for e in edges:
        s = e.get('source')
        t = e.get('target')
        r = e.get('relation')
        relations[r] += 1
        edge_keys[(s, t, r)] += 1
        if r in illegal_relations:
            illegal.append(e)
        desc = str(e.get('description') or e.get('desc') or '').strip()
        if len(desc) < 10:
            no_desc.append(e)
        if s == t:
            selfloops.append(e)
        if s not in node_ids or t not in node_ids:
            orphan.append({
                'source': s,
                'target': t,
                'relation': r,
                'source_exists': s in node_ids,
                'target_exists': t in node_ids,
                'description': desc[:120],
            })
    dup_edges = [
        {'source': k[0], 'target': k[1], 'relation': k[2], 'count': v}
        for k, v in edge_keys.items() if v > 1
    ]
    REPORT['edges'] = {
        'rows': len(edges),
        'relations': dict(relations),
        'illegal_nonbio_count': len(illegal),
        'missing_description_count': len(no_desc),
        'orphan_count': len(orphan),
        'selfloops_count': len(selfloops),
        'duplicate_edge_count': len(dup_edges),
        'illegal_examples': illegal[:10],
        'no_desc_examples': no_desc[:10],
        'orphan_examples': orphan[:20],
        'selfloop_examples': selfloops[:10],
        'duplicate_examples': dup_edges[:20],
    }
    if illegal:
        add_issue('error', 'Illegal non-biological edge relations remain', illegal[:5])
    if no_desc:
        add_issue('warning', 'Edges missing useful description', no_desc[:5])
    if orphan:
        add_issue('warning', 'Edges reference missing nodes', orphan[:10])
    if selfloops:
        add_issue('warning', 'Self-loop edges exist', selfloops[:5])
    if dup_edges:
        add_issue('warning', 'Duplicate edges exist', dup_edges[:10])


def audit_fulltext_cache():
    fc = ROOT / 'fulltext_cache'
    if not fc.exists() and (ROOT / 'fulltext').exists():
        fc = ROOT / 'fulltext'
    if not fc.exists():
        add_issue('warning', 'fulltext/fulltext_cache missing')
        return
    txts = sorted(fc.glob('*.txt'))
    metas = sorted(fc.glob('*.metadata.json'))
    txt_stems = {p.name[:-4] for p in txts}
    meta_stems = {p.name[:-14] for p in metas}
    naming = []
    tiny = []
    cloudflare = []
    bad_meta = []
    txt_no_meta = []
    meta_no_txt = []
    for p in txts:
        stem = p.name[:-4]
        if not re.match(r'^(PMID_\d+|PMC\d+|ARXIV_[A-Za-z0-9_.-]+|BIORXIV_10\.[A-Za-z0-9_.-]+|MEDRXIV_10\.[A-Za-z0-9_.-]+)', stem):
            naming.append(p.name)
        sz = p.stat().st_size
        text = p.read_text(encoding='utf-8-sig', errors='replace')
        if sz < 1000:
            tiny.append({'file': p.name, 'size': sz, 'preview': text[:200]})
        if re.search(r'Just a moment|Verify you are human|Cloudflare', text[:2000], re.I) or ('security verification' in text[:2000].lower()):
            cloudflare.append({'file': p.name, 'preview': text[:300]})
        if stem.startswith(('BIORXIV_', 'MEDRXIV_')) and stem not in meta_stems:
            txt_no_meta.append(p.name)
    for p in metas:
        stem = p.name[:-14]
        try:
            m = json.loads(p.read_text(encoding='utf-8-sig'))
            doi = str(m.get('doi') or m.get('extractedDoi') or (m.get('apiRecord') or {}).get('doi') or '')
            title = str(m.get('title') or '')
            if not (doi or title):
                bad_meta.append({'file': p.name, 'reason': 'no doi/title'})
        except Exception as e:
            bad_meta.append({'file': p.name, 'reason': repr(e)})
        if stem not in txt_stems:
            meta_no_txt.append(p.name)
    hashes = defaultdict(list)
    for p in txts:
        hashes[hashlib.sha256(p.read_bytes()).hexdigest()].append(p.name)
    duplicate_content = [v for v in hashes.values() if len(v) > 1]
    total_size = sum(p.stat().st_size for p in fc.rglob('*') if p.is_file())
    REPORT['fulltext_cache'] = {
        'txt_files': len(txts),
        'metadata_files': len(metas),
        'dirs': [p.name for p in fc.iterdir() if p.is_dir()],
        'total_size': fmt(total_size),
        'naming_issues': naming,
        'tiny_files': tiny,
        'cloudflare_residue': cloudflare,
        'bio_preprint_txt_without_meta': txt_no_meta,
        'metadata_without_txt': meta_no_txt,
        'bad_metadata': bad_meta,
        'duplicate_content': duplicate_content,
    }
    if naming:
        add_issue('warning', 'Fulltext files with nonstandard names', naming[:10])
    if tiny:
        add_issue('warning', 'Tiny fulltext cache files may be abstracts only', tiny[:10])
    if cloudflare:
        add_issue('error', 'Cloudflare/security pages remain in fulltext cache', cloudflare)
    if bad_meta:
        add_issue('warning', 'Bad metadata files', bad_meta[:10])
    if duplicate_content:
        add_issue('warning', 'Duplicate fulltext contents', duplicate_content[:10])


def write_report():
    today = date.today().isoformat()
    out_json = ROOT / 'daily_digest' / f'selfcheck_{today}.json'
    out_md = ROOT / 'daily_digest' / f'selfcheck_{today}.md'
    out_json.parent.mkdir(exist_ok=True)
    out_json.write_text(json.dumps(REPORT, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
    lines = []
    lines.append(f'# Literature knowledge graph self-check {today}')
    lines.append('')
    lines.append(f"Timestamp: {REPORT['timestamp']}")
    lines.append(f"Total: {REPORT['inventory']['total_files']} files, {REPORT['inventory']['total_size_human']}")
    lines.append('')
    lines.append('## Inventory')
    for x in REPORT['inventory']['top_level']:
        lines.append(f"- {x['name']}: {x['files']} files, {x['size']}")
    lines.append('')
    lines.append('## DB integrity')
    for name in ['papers.db', 'edges.db', 'concepts.db']:
        d = REPORT['db'][name]
        lines.append(f"- {name}: rows={d['rows']}, bad_lines={d['bad_lines']}, size={d['size']}")
    pi = REPORT['db']['papers_integrity']
    lines.append(f"- papers unique IDs: {pi['unique_ids']}; duplicate IDs: {len(pi['duplicate_ids'])}; missing IDs: {pi['missing_ids']}")
    lines.append('')
    lines.append('## S-tier')
    st = REPORT['s_tier']
    lines.append(f"- total_S={st['total_S']}; empty_s={st['empty_s_count']}; weak_s={st['weak_s_count']}")
    lines.append(f"- field pass counts: {st['field_pass_counts']}")
    lines.append('')
    lines.append('## Edges')
    ed = REPORT['edges']
    lines.append(f"- rows={ed['rows']}; illegal_nonbio={ed['illegal_nonbio_count']}; missing_desc={ed['missing_description_count']}; orphan={ed['orphan_count']}; selfloops={ed['selfloops_count']}; duplicates={ed['duplicate_edge_count']}")
    lines.append('')
    lines.append('## Fulltext cache')
    ft = REPORT['fulltext_cache']
    lines.append(f"- txt={ft.get('txt_files')}; metadata={ft.get('metadata_files')}; size={ft.get('total_size')}; dirs={ft.get('dirs')}")
    lines.append(f"- cloudflare_residue={len(ft.get('cloudflare_residue', []))}; tiny_files={len(ft.get('tiny_files', []))}; naming_issues={len(ft.get('naming_issues', []))}; bad_metadata={len(ft.get('bad_metadata', []))}")
    lines.append('')
    lines.append('## Issues')
    for x in REPORT['issues']:
        lines.append(f"- [{x['severity']}] {x['message']}: {str(x.get('data', ''))[:500]}")
    lines.append('')
    lines.append('## Warnings')
    for x in REPORT['warnings'][:80]:
        lines.append(f"- [{x['severity']}] {x['message']}: {str(x.get('data', ''))[:500]}")
    out_md.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return out_json, out_md


def main():
    scan_inventory()
    check_forbidden_and_ports()
    papers = load_db('papers.db')
    edges = load_db('edges.db')
    concepts = load_db('concepts.db')
    audit_papers(papers)
    audit_s_tier(papers)
    audit_concepts(concepts)
    audit_edges(edges, papers, concepts)
    audit_fulltext_cache()
    out_json, out_md = write_report()
    print(json.dumps({
        'report_json': str(out_json),
        'report_md': str(out_md),
        'issues': len(REPORT['issues']),
        'warnings': len(REPORT['warnings']),
        'summary': {
            'files': REPORT['inventory']['total_files'],
            'size': REPORT['inventory']['total_size_human'],
            'db_bad_lines': {k: REPORT['db'][k]['bad_lines'] for k in ['papers.db', 'edges.db', 'concepts.db']},
            'paper_duplicate_ids': len(REPORT['db']['papers_integrity']['duplicate_ids']),
            's_empty': REPORT['s_tier']['empty_s_count'],
            's_weak': REPORT['s_tier']['weak_s_count'],
            'edge_illegal': REPORT['edges']['illegal_nonbio_count'],
            'edge_no_desc': REPORT['edges']['missing_description_count'],
            'edge_orphan': REPORT['edges']['orphan_count'],
            'edge_selfloops': REPORT['edges']['selfloops_count'],
            'fulltext_cloudflare': len(REPORT['fulltext_cache'].get('cloudflare_residue', [])),
            'fulltext_tiny': len(REPORT['fulltext_cache'].get('tiny_files', [])),
        }
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
