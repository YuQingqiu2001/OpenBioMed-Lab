#!/usr/bin/env python3
"""Build a portable, self-contained HTML view of the literature graph."""

import json
from pathlib import Path
from typing import Any

from workspace_paths import KG_ROOT, ensure_workspace


PAPERS_DB = KG_ROOT / "papers.db"
CONCEPTS_DB = KG_ROOT / "concepts.db"
EDGES_DB = KG_ROOT / "edges.db"
OUTPUT = KG_ROOT / "network.html"

GROUP_COLORS = {
    "pubmed": "#4ECDC4",
    "arxiv": "#FF6B6B",
    "biorxiv": "#45B7D1",
    "medrxiv": "#96CEB4",
    "concept": "#FFD93D",
    "unknown": "#8B949E",
}


def load_ndjson(path: Path) -> list[dict[str, Any]]:
    records = []
    if not path.exists():
        return records
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                print(f"Skipping invalid JSON in {path.name}:{line_number}: {exc}")
                continue
            if isinstance(record, dict):
                records.append(record)
    return records


def detect_group(identifier: str, source: str = "") -> str:
    value = f"{identifier} {source}".lower()
    if "pubmed" in value or identifier.upper().startswith("PMID:"):
        return "pubmed"
    if "biorxiv" in value:
        return "biorxiv"
    if "medrxiv" in value:
        return "medrxiv"
    if "arxiv" in value or identifier.upper().startswith("ARXIV:"):
        return "arxiv"
    if "concept" in value or identifier.upper().startswith("CONCEPT:"):
        return "concept"
    return "unknown"


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def text_value(value: Any, limit: int = 180) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        text = value
    elif isinstance(value, dict):
        preferred = (
            value.get("summary")
            or value.get("claim")
            or value.get("finding")
            or value.get("description")
        )
        text = str(preferred) if preferred else json.dumps(value, ensure_ascii=False)
    elif isinstance(value, (list, tuple)):
        text = "; ".join(text_value(item, limit=limit) for item in value)
    else:
        text = str(value)
    return " ".join(text.split())[:limit]


def paper_analysis(paper: dict[str, Any]) -> dict[str, Any]:
    analysis = paper.get("analysis")
    return analysis if isinstance(analysis, dict) else {}


def paper_claims(paper: dict[str, Any]) -> list[Any]:
    analysis = paper_analysis(paper)
    tier2 = analysis.get("tier2_full")
    tier2 = tier2 if isinstance(tier2, dict) else {}
    return as_list(
        paper.get("claims")
        or paper.get("tier3_ces_chains")
        or analysis.get("ces_chains")
        or tier2.get("key_findings")
    )


def paper_subquestions(paper: dict[str, Any]) -> list[Any]:
    analysis = paper_analysis(paper)
    return as_list(
        paper.get("tier2_subquestions")
        or analysis.get("subquestions")
    )


def paper_finding(paper: dict[str, Any], claims: list[Any]) -> str:
    analysis = paper_analysis(paper)
    candidates = (
        paper.get("core_findings"),
        paper.get("tier2_core_question"),
        analysis.get("core_question"),
        claims,
        paper.get("abstract"),
    )
    for candidate in candidates:
        text = text_value(candidate)
        if text:
            return text
    return ""


def build_graph() -> dict[str, list[dict[str, Any]]]:
    papers = load_ndjson(PAPERS_DB)
    concepts = load_ndjson(CONCEPTS_DB)
    edges = load_ndjson(EDGES_DB)
    nodes: list[dict[str, Any]] = []
    links: list[dict[str, Any]] = []
    seen: set[str] = set()

    for paper in papers:
        identifier = str(paper.get("id") or paper.get("paper_id") or "")
        if not identifier or identifier in seen:
            continue
        seen.add(identifier)
        claims = paper_claims(paper)
        subquestions = paper_subquestions(paper)
        analysis = paper_analysis(paper)
        explicit_tier = paper.get("analysis_tier") or analysis.get("tier")
        tier = str(explicit_tier or ("S" if claims else "?"))
        deep = len(claims) >= 3 and len(subquestions) >= 3
        group = detect_group(identifier, str(paper.get("source") or ""))
        size = 18 if tier.upper() == "S" and deep else 12 if tier.upper() == "S" else 7
        nodes.append(
            {
                "id": identifier,
                "label": text_value(paper.get("title") or "Untitled", limit=60),
                "group": group,
                "color": GROUP_COLORS[group],
                "size": size,
                "claims": len(claims),
                "deep": deep,
                "tier": tier,
                "journal": text_value(paper.get("journal"), limit=80),
                "year": text_value(paper.get("year"), limit=8),
                "findings": paper_finding(paper, claims),
            }
        )

    for concept in concepts:
        identifier = str(concept.get("id") or concept.get("concept_id") or "")
        if not identifier or identifier in seen:
            continue
        seen.add(identifier)
        nodes.append(
            {
                "id": identifier,
                "label": text_value(concept.get("name") or identifier, limit=40),
                "group": "concept",
                "color": GROUP_COLORS["concept"],
                "size": 14,
                "claims": 0,
                "deep": True,
                "tier": "CONCEPT",
                "journal": text_value(concept.get("type"), limit=80),
                "year": "",
                "findings": text_value(concept.get("definition")),
            }
        )

    for edge in edges:
        source = str(edge.get("source") or edge.get("from") or "")
        target = str(edge.get("target") or edge.get("to") or "")
        if not source or not target:
            continue
        relation = text_value(edge.get("relation") or edge.get("type") or "related", 80)
        links.append(
            {
                "source": source,
                "target": target,
                "relation": relation,
                "description": text_value(edge.get("description"), 120),
                "isConcept": detect_group(source) == "concept"
                or detect_group(target) == "concept",
            }
        )

    return {"nodes": nodes, "links": links}


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Biomedical Literature Network</title>
<style>
*{box-sizing:border-box} body{margin:0;overflow:hidden;background:#0d1117;color:#c9d1d9;font:13px system-ui,sans-serif}
canvas{display:block;width:100vw;height:100vh}.panel{position:fixed;z-index:2;padding:12px 14px;background:#161b22e8;border:1px solid #30363d;border-radius:8px}
#legend{top:12px;left:12px}#stats{top:12px;right:12px}.row{margin:4px 0}.dot{display:inline-block;width:10px;height:10px;margin-right:7px;border-radius:50%}
#tip{display:none;left:12px;bottom:12px;max-width:480px;pointer-events:none}.muted{color:#8b949e}.title{color:#58a6ff;font-weight:700}
</style>
</head>
<body>
<canvas id="graph"></canvas>
<div id="legend" class="panel">
  <div class="title">Sources</div>
  <div class="row"><span class="dot" style="background:#4ECDC4"></span>PubMed</div>
  <div class="row"><span class="dot" style="background:#FF6B6B"></span>arXiv</div>
  <div class="row"><span class="dot" style="background:#45B7D1"></span>bioRxiv</div>
  <div class="row"><span class="dot" style="background:#96CEB4"></span>medRxiv</div>
  <div class="row"><span class="dot" style="background:#FFD93D"></span>Concept</div>
</div>
<div id="stats" class="panel">
  <div class="title">Knowledge Graph</div>
  <div class="row">Papers: <b>__PAPERS__</b></div>
  <div class="row">Concepts: <b>__CONCEPTS__</b></div>
  <div class="row">S-tier: <b>__S_TIER__</b></div>
  <div class="row">Edges: <b>__EDGES__</b></div>
  <div class="row muted">Drag to pan, wheel to zoom</div>
</div>
<div id="tip" class="panel"><div class="title"></div><div class="muted"></div><div class="body"></div></div>
<script>
const DATA={nodes:__NODES__,links:__LINKS__};
const canvas=document.getElementById("graph"),ctx=canvas.getContext("2d"),tip=document.getElementById("tip");
let width=0,height=0,dpr=1,scale=1,offsetX=0,offsetY=0,dragging=false,lastX=0,lastY=0;
const nodes=DATA.nodes,links=DATA.links,nodeMap=new Map(nodes.map(n=>[n.id,n]));
function resize(){dpr=window.devicePixelRatio||1;width=innerWidth;height=innerHeight;canvas.width=width*dpr;canvas.height=height*dpr;canvas.style.width=width+"px";canvas.style.height=height+"px"}
function seed(){nodes.forEach((n,i)=>{const a=i/Math.max(nodes.length,1)*Math.PI*2,r=Math.min(width,height)*(.15+.25*(i%5)/5);n.x=width/2+Math.cos(a)*r;n.y=height/2+Math.sin(a)*r;n.vx=0;n.vy=0})}
function simulate(){for(let i=0;i<nodes.length;i++){for(let j=i+1;j<nodes.length;j++){const a=nodes[i],b=nodes[j],dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(dx*dx+dy*dy,100),f=900/d2;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f}}
for(const link of links){const a=nodeMap.get(link.source),b=nodeMap.get(link.target);if(!a||!b)continue;const dx=b.x-a.x,dy=b.y-a.y,d=Math.max(Math.hypot(dx,dy),1),f=(d-110)*.0015;a.vx+=dx/d*f;a.vy+=dy/d*f;b.vx-=dx/d*f;b.vy-=dy/d*f}
for(const n of nodes){n.vx+=(width/2-n.x)*.0002;n.vy+=(height/2-n.y)*.0002;n.vx*=.86;n.vy*=.86;n.x+=n.vx;n.y+=n.vy}}
function draw(){ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);ctx.save();ctx.translate(offsetX,offsetY);ctx.scale(scale,scale);
ctx.lineWidth=1;for(const link of links){const a=nodeMap.get(link.source),b=nodeMap.get(link.target);if(!a||!b)continue;ctx.strokeStyle=link.isConcept?"#FFD93Daa":"#58a6ff66";ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke()}
ctx.font="11px system-ui";for(const n of nodes){ctx.fillStyle=n.color;ctx.strokeStyle=n.deep?"#fff":n.color;ctx.lineWidth=n.deep?2:1;ctx.beginPath();ctx.arc(n.x,n.y,n.size,0,Math.PI*2);ctx.fill();ctx.stroke();if(n.size>=12){ctx.fillStyle="#c9d1d9";ctx.fillText(n.label.slice(0,28),n.x+n.size+4,n.y+4)}}ctx.restore()}
function loop(){simulate();draw();requestAnimationFrame(loop)}
function screenToGraph(x,y){return{x:(x-offsetX)/scale,y:(y-offsetY)/scale}}
canvas.addEventListener("wheel",e=>{e.preventDefault();const before=screenToGraph(e.clientX,e.clientY),factor=e.deltaY>0?.9:1.1;scale=Math.max(.2,Math.min(5,scale*factor));offsetX=e.clientX-before.x*scale;offsetY=e.clientY-before.y*scale},{passive:false});
canvas.addEventListener("mousedown",e=>{dragging=true;lastX=e.clientX;lastY=e.clientY});
addEventListener("mouseup",()=>dragging=false);addEventListener("mousemove",e=>{if(dragging){offsetX+=e.clientX-lastX;offsetY+=e.clientY-lastY;lastX=e.clientX;lastY=e.clientY;tip.style.display="none";return}
const p=screenToGraph(e.clientX,e.clientY);let hit=null;for(const n of nodes){if(Math.hypot(n.x-p.x,n.y-p.y)<=n.size+4){hit=n;break}}if(!hit){tip.style.display="none";return}tip.style.display="block";tip.querySelector(".title").textContent=(hit.tier==="S"?"[S] ":"")+hit.label;tip.querySelector(".muted").textContent=[hit.group,hit.journal,hit.year,"claims: "+hit.claims].filter(Boolean).join(" | ");tip.querySelector(".body").textContent=hit.findings||""});
addEventListener("resize",()=>{resize()});resize();seed();loop();
</script>
</body>
</html>
"""


def render(graph: dict[str, list[dict[str, Any]]]) -> str:
    nodes = graph["nodes"]
    links = graph["links"]
    nodes_json = json.dumps(nodes, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    links_json = json.dumps(links, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    values = {
        "__PAPERS__": str(sum(node["group"] != "concept" for node in nodes)),
        "__CONCEPTS__": str(sum(node["group"] == "concept" for node in nodes)),
        "__S_TIER__": str(sum(str(node.get("tier", "")).upper() == "S" for node in nodes)),
        "__EDGES__": str(len(links)),
        "__NODES__": nodes_json,
        "__LINKS__": links_json,
    }
    html = HTML_TEMPLATE
    for marker, value in values.items():
        html = html.replace(marker, value)
    return html


def main() -> None:
    ensure_workspace()
    graph = build_graph()
    OUTPUT.write_text(render(graph), encoding="utf-8")
    print(
        f"Wrote {OUTPUT} with {len(graph['nodes'])} nodes "
        f"and {len(graph['links'])} edges."
    )


if __name__ == "__main__":
    main()
