"""Freeze the arithmetic mean of all same-cancer, same-reference teachers.

Deployment targets, deliberately NOT out-of-fold biological validation.
Reads only H&E feature caches, checkpoint states and geometry.
"""
from pathlib import Path
import argparse, json, sys, time, hashlib, os
import numpy as np, pandas as pd, torch
import os
from pathlib import Path
ROOT=Path(os.environ['COVARST_STUDY_ROOT']).resolve()
OFF=ROOT
OUT=Path(os.environ['COVARST_BUILD_ROOT']).resolve()

sys.path.insert(0,str(Path(__file__).parent))
from distill_final_models import construct, dump, sha, OUT, _to_device_graph_pair
from run_virchow2_fseg_k64_mapper import BoundedINR
from openst_final.he_uni_mamba_bnn_inr import normalize_slide_coordinates

def main():
 p=argparse.ArgumentParser();p.add_argument('--tissue',action='append');a=p.parse_args()
 torch.set_num_threads(8);torch.set_float32_matmul_precision('high')
 torch.cuda.set_per_process_memory_fraction(6*(1<<30)/torch.cuda.get_device_properties(0).total_memory)
 dest=OUT/'ensemble_targets_v1';dest.mkdir(exist_ok=True)
 (dest/'pid.txt').write_text(str(os.getpid()))
 registry=pd.read_csv(OFF/'checkpoints/checkpoint_registry.csv')
 tissues=a.tissue or sorted(registry.tissue_key.unique())
 for tissue in tissues:
  folder=dest/tissue;folder.mkdir(exist_ok=True)
  rows=json.loads((OUT/'teacher_cache'/tissue/'manifest.json').read_text())
  rr=registry.loc[registry.tissue_key==tissue].sort_values('fold_index')
  contract={'method':'equal_weight_probability_ensemble_v1','teacher_folds':len(rr),'tissue':tissue,'reference_sha256':sha(OFF/'references'/tissue/'consensus_reference_probability_10k.npy'),'gene_order_sha256':sha(OFF/'references'/tissue/'gene_order_10k.npy'),'teacher_checkpoint_sha256':rr.checkpoint_sha256.tolist(),'input':'H&E features and coordinates only','aggregation':'arithmetic mean of softmax Program probabilities after accepted INR correction','student_validation':'deployment teacher fidelity only; not independent biological validation','code_sha256':sha(__file__)}
  if (folder/'contract.json').exists():assert json.loads((folder/'contract.json').read_text())==contract
  else:dump(folder/'contract.json',contract)
  if (folder/'complete.json').exists():print(json.dumps({'reused':tissue}),flush=True);continue
  started=time.time();teachers=[]
  for r in rr.itertuples():
   path=OFF/str(r.checkpoint_relative_path).replace('\\','/')
   assert sha(path)==r.checkpoint_sha256
   cp=torch.load(path,map_location='cpu',weights_only=False)
   if cp.get('consensus_reference_sha256'):assert cp['consensus_reference_sha256']==contract['reference_sha256']
   groups=np.load(path.parent/'program_group_index.npy',allow_pickle=False)
   model=construct(cp,groups,'cuda').eval();inr=None
   if cp.get('inr_accepted') and cp.get('inr') is not None:
    inr=BoundedINR(int(cp['hidden']),int(cp['programs'])*2).cuda().eval();inr.load_state_dict(cp['inr'],strict=True)
   teachers.append((model,inr,float(cp.get('support_gamma',1.)),int(cp['programs'])))
  results=[]
  with torch.inference_mode():
   for j,row in enumerate(rows):
    target=folder/(Path(row['cache']).stem+'.npz')
    if not target.exists():
     with np.load(row['cache'],allow_pickle=False) as z:
      f=torch.as_tensor(z['features'].astype(np.float32),device='cuda')
      gs=_to_device_graph_pair(tuple(tuple(z[f'{s}_{k}'] for k in ('source','target','prior')) for s in ('local','context')),'cuda')
      coord=torch.as_tensor(normalize_slide_coordinates(z['coords']),device='cuda')
     sum_r=torch.zeros((len(f),row['programs']),device='cuda');sum_c=sum_r.clone();square_r=sum_r.clone()
     for model,inr,gamma,k in teachers:
      r,c,context,_=model(f,gs,support_gamma=gamma,domain_strength=0.)
      if inr is not None:
       dr,dc=inr(coord,context).split(k,dim=1);r=r+dr;c=c+dc
      r=torch.softmax(r,1);c=torch.softmax(c,1);sum_r+=r;sum_c+=c;square_r+=r.square()
     mean_r=sum_r/len(teachers);mean_c=sum_c/len(teachers)
     assert torch.isfinite(mean_r).all() and torch.allclose(mean_r.sum(1),torch.ones(len(f),device='cuda'),atol=2e-6)
     with target.with_suffix('.partial').open('wb') as h:
      np.savez(h,teacher_rna=mean_r.cpu().numpy(),teacher_composition=mean_c.cpu().numpy(),disagreement_mean_square=float((square_r/len(teachers)-mean_r.square()).sum(1).mean()))
     target.with_suffix('.partial').replace(target)
    record={**row,'target':str(target),'target_sha256':sha(target),'teacher_fold':'all','teacher_ensemble_folds':len(teachers),'teacher_patient_held_out':False}
    results.append(record)
    if (j+1)%20==0 or j+1==len(rows):print(json.dumps({'tissue':tissue,'slides':j+1,'total':len(rows),'seconds':round(time.time()-started,1)}),flush=True)
  dump(folder/'manifest.json',results);dump(folder/'complete.json',{'tissue':tissue,'slides':len(rows),'spots':sum(r['spots'] for r in rows),'teachers':len(teachers),'seconds':time.time()-started})
  del teachers;torch.cuda.empty_cache()

if __name__=='__main__':main()
