"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from pathlib import Path

import argparse, copy, hashlib, json, math, os, sys, time

import numpy as np, pandas as pd, torch

from train_virchow2_pan_atlas_consensus_mapper import PanAtlasDualMapper, _to_device_graph_pair

EPS=1e-8

def sha(p):
 d=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):d.update(b)
 return d.hexdigest()

def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.partial');q.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False));q.replace(p)

def save(p,x):
 q=p.with_suffix('.partial');torch.save(x,q);q.replace(p)

def js(p,q):
 m=(p+q)*.5
 return .5*((p*(p.clamp_min(EPS).log()-m.clamp_min(EPS).log())).sum(1)+(q*(q.clamp_min(EPS).log()-m.clamp_min(EPS).log())).sum(1)).mean()

def corr(a,b):
 a=a-a.mean(0);b=b-b.mean(0);den=(a.square().sum(0)*b.square().sum(0)).sqrt()
 ok=den>1e-8
 return 1-(a[:,ok]*b[:,ok]).sum(0).div(den[ok]).mean() if ok.any() else a.sum()*0

def construct(cp,groups,device):
 m=PanAtlasDualMapper(int(cp['input_dim']),hidden=int(cp['hidden']),programs=int(cp['programs']),program_groups=groups,graph_layers=int(cp['graph_layers']),context_graph_layers=int(cp['context_graph_layers']),graph_heads=int(cp['graph_heads']),dropout=.10,type_lift_alpha=float(cp['type_lift_alpha']),sources=int(cp['sources']),technologies=int(cp['technologies']),support_mass=float(cp['support_mass']),hierarchical_logits=bool(cp.get('hierarchical_logits',False)))
 m.load_state_dict(cp['model'],strict=True);return m.to(device)

def load(row,device):
 with np.load(row['cache'],allow_pickle=False) as z:
  f=torch.as_tensor(z['features'].astype(np.float32),device=device)
  graph=tuple(tuple(z[f'{prefix}_{k}'] for k in ('source','target','prior')) for prefix in ('local','context'))
  graphs=_to_device_graph_pair(graph,device)
  r=torch.as_tensor(z['teacher_rna'],device=device);c=torch.as_tensor(z['teacher_composition'],device=device)
 return f,graphs,r,c

def basis(row,reference,bybatch,batches,device):
 key=row['source']+'|'+row['technology']
 if key in batches:return torch.as_tensor(np.array(bybatch[batches[key]],np.float32),device=device)
 return torch.as_tensor(np.array(reference,np.float32),device=device)

def probabilities(model,features,graphs,cp):
 r,c,_,aux=model(features,graphs,support_gamma=float(cp.get('support_gamma',1.0)),domain_strength=0.0)
 return torch.softmax(r,1),torch.softmax(c,1)

@torch.inference_mode()
def evaluate(model,rows,cp,ref,bybatch,batches,device,detail=False):
 model.eval();records=[]
 for row in rows:
  f,graphs,tr,tc=load(row,device);pr,pc=probabilities(model,f,graphs,cp);w=basis(row,ref,bybatch,batches,device)
  values=[];huber=0.;dot=torch.zeros(len(f),device=device);n1=dot.clone();n2=dot.clone()
  for start in range(0,w.shape[1],512):
   b=w[:,start:start+512];pred=torch.log1p(10000*(pr@b));target=torch.log1p(10000*(tr@b))
   x=pred-pred.mean(0);y=target-target.mean(0);den=(x.square().sum(0)*y.square().sum(0)).sqrt();ok=den>1e-8
   values.extend(((x[:,ok]*y[:,ok]).sum(0)/den[ok]).cpu().numpy().tolist())
   huber+=float(torch.nn.functional.smooth_l1_loss(pred,target))*(b.shape[1]/w.shape[1])
   dot+=(pred*target).sum(1);n1+=pred.square().sum(1);n2+=target.square().sum(1)
  pcc=float(np.median(values)) if values else None
  src,tgt,_=graphs[0];d1=pr[src]-pr[tgt];d2=tr[src]-tr[tgt]
  edge_den=float(torch.linalg.vector_norm(d1)*torch.linalg.vector_norm(d2));edge_cos=float((d1*d2).sum())/edge_den if edge_den>1e-10 else None
  rjs=float(js(pr,tr));cjs=float(js(pc,tc))
  record={k:row[k] for k in ('tissue','slide','patient','source','spots','teacher_fold')}
  record.update(rna_jsd=rjs,composition_jsd=cjs,expression_gene_pcc_median=pcc,finite_expression_genes=len(values),expression_log_huber=huber,expression_spot_cosine_median=float(torch.median(dot/(n1*n2).sqrt().clamp_min(EPS))),spatial_theta_gradient_cosine=edge_cos,selection_loss=.7*rjs+.45*cjs+2*huber+4*(1-pcc if pcc is not None else 1))
  records.append(record)
 frame=pd.DataFrame(records)
 # Equal patients inside each source, then equal source weights.
 balanced=frame.groupby(['source','patient']).mean(numeric_only=True).groupby('source').mean(numeric_only=True).mean(numeric_only=True)
 summary={k:float(balanced[k]) for k in ['rna_jsd','composition_jsd','expression_gene_pcc_median','expression_log_huber','expression_spot_cosine_median','spatial_theta_gradient_cosine','selection_loss'] if np.isfinite(balanced[k])}
 summary.update(slides=len(rows),patients=len({r['patient'] for r in rows}),spots=sum(r['spots'] for r in rows))
 return summary,frame

def seed(n):
 np.random.seed(n);torch.manual_seed(n)
 if torch.cuda.is_available():torch.cuda.manual_seed_all(n)


import os
from pathlib import Path
ROOT=Path(os.environ['COVARST_STUDY_ROOT']).resolve()
OFF=ROOT
OUT=Path(os.environ['COVARST_BUILD_ROOT']).resolve()
