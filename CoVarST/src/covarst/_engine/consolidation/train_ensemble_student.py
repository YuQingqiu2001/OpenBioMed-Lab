"""Fit one unchanged-capacity student per cancer to frozen ensemble targets.

All eligible slides enter deployment consolidation. Metrics here are teacher
fidelity, never a replacement for the original patient-out biological scores.
"""
from pathlib import Path
import argparse, copy, json, os, time
import numpy as np,pandas as pd,torch
import distill_final_models as base
ROOT=base.ROOT;OFF=base.OFF;OUT=base.OUT
CPU={}

def load(row,device):
 key=row['cache']
 if key not in CPU:
  with np.load(key,allow_pickle=False) as z:
   f=z['features'].astype(np.float32)
   graph=tuple(tuple(z[f'{prefix}_{k}'] for k in ('source','target','prior')) for prefix in ('local','context'))
  with np.load(row['target'],allow_pickle=False) as z:r=z['teacher_rna'];c=z['teacher_composition']
  CPU[key]=(f,graph,r,c)
 f,g,r,c=CPU[key]
 return torch.as_tensor(f,device=device),base._to_device_graph_pair(g,device),torch.as_tensor(r,device=device),torch.as_tensor(c,device=device)
base.load=load

def step(model,opt,row,cp,device):
 f,g,tr,tc=load(row,device);opt.zero_grad(set_to_none=True)
 pr,pc=base.probabilities(model,f,g,cp)
 # Match probabilities, spatial variation and per-Program centred signals.
 src,tgt,_=g[0];dr=pr[src]-pr[tgt];dt=tr[src]-tr[tgt]
 scale=tr.shape[1]
 loss=2*base.js(pr,tr)+base.js(pc,tc)+scale*torch.nn.functional.mse_loss(pr,tr)+scale*torch.nn.functional.mse_loss(pc,tc)+base.corr(pr,tr)
 loss+=.5*scale*torch.nn.functional.mse_loss(dr,dt)
 loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5);opt.step()
 return float(loss.detach())

def run(tissue,a):
 targetroot=OUT/'ensemble_targets_v1'/tissue
 if not (targetroot/'complete.json').exists():return False
 rows=json.loads((targetroot/'manifest.json').read_text());dest=a.output/tissue;dest.mkdir(parents=True,exist_ok=True)
 if (dest/'completed.json').exists():return True
 rr=pd.read_csv(OFF/'checkpoints/checkpoint_registry.csv');rr=rr[rr.tissue_key==tissue].sort_values('fold_index')
 entry=rr.iloc[0];path=OFF/str(entry.checkpoint_relative_path).replace('\\','/')
 cp=torch.load(path,map_location='cpu',weights_only=False);groups=np.load(path.parent/'program_group_index.npy',allow_pickle=False)
 refroot=OFF/'references'/tissue;ref=np.load(refroot/'consensus_reference_probability_10k.npy',mmap_mode='r');bybatch=np.load(refroot/'consensus_reference_probability_by_batch_10k.npy',mmap_mode='r');batches={str(r.batch):int(r.batch_index) for r in pd.read_csv(refroot/'batch_adapter_metadata.csv').itertuples()}
 contract={'method':'same_capacity_ensemble_probability_distillation_v1','tissue':tissue,'epochs_max':a.epochs,'learning_rate':a.lr,'seed':a.seed,'teacher_folds':len(rr),'target_contract':json.loads((targetroot/'contract.json').read_text()),'reference_sha256':base.sha(refroot/'consensus_reference_probability_10k.npy'),'gene_order_sha256':base.sha(refroot/'gene_order_10k.npy'),'trainer_sha256':base.sha(__file__),'architecture_unchanged':True,'parameter_count_unchanged':True,'initial_teacher_sha256':base.sha(path),'all_eligible_spots_used':sum(r['spots'] for r in rows),'validation_interpretation':'same-cancer deployment ensemble imitation on consolidation inputs, not independent biological or clinical validation','teacher_ensemble_may_have_seen_consolidation_patients':True,'student_raw_RNA_read':False,'checkpoint_scope':'H&E-to-Program mapper only; fixed reference matrix distributed separately','student_INR':False,'gate':{'rna_jsd_max':.05,'composition_jsd_max':.05,'expression_gene_pcc_min':.90,'expression_spot_cosine_min':.98,'spatial_gradient_cosine_min':.85}}
 if (dest/'contract.json').exists():assert json.loads((dest/'contract.json').read_text())==contract
 else:base.dump(dest/'contract.json',contract)
 device=torch.device('cuda');base.seed(a.seed);model=base.construct(cp,groups,device)
 # Disabled stochastic dropout during imitation; modules/capacity stay identical.
 for m in model.modules():
  if isinstance(m,torch.nn.Dropout):m.p=0.
 opt=torch.optim.AdamW(model.parameters(),lr=a.lr,weight_decay=1e-4)
 history=[];start=1;best=float('inf');rng=np.random.default_rng(a.seed)
 if (dest/'resume.pt').exists():
  s=torch.load(dest/'resume.pt',map_location=device,weights_only=False);model.load_state_dict(s['model']);opt.load_state_dict(s['optimizer']);history=s['history'];start=s['epoch']+1;best=s['best'];rng.bit_generator.state=s['rng']
 if not (dest/'initial_fidelity.json').exists():
  metric,_=base.evaluate(model,rows,cp,ref,bybatch,batches,device);base.dump(dest/'initial_fidelity.json',metric)
 for epoch in range(start,a.epochs+1):
  t=time.time();model.eval();losses=[]
  for index in rng.permutation(len(rows)):losses.append(step(model,opt,rows[int(index)],cp,device))
  rec={'epoch':epoch,'loss':float(np.mean(losses)),'seconds':time.time()-t}
  passed=False
  if epoch%5==0 or epoch==a.epochs:
   metric,frame=base.evaluate(model,rows,cp,ref,bybatch,batches,device)
   rec.update(metric);gate=contract['gate']
   passed=metric['rna_jsd']<=gate['rna_jsd_max'] and metric['composition_jsd']<=gate['composition_jsd_max'] and metric['expression_gene_pcc_median']>=gate['expression_gene_pcc_min'] and metric['expression_spot_cosine_median']>=gate['expression_spot_cosine_min'] and metric.get('spatial_theta_gradient_cosine',-1)>=gate['spatial_gradient_cosine_min']
   if metric['selection_loss']<best:
    best=metric['selection_loss'];base.save(dest/'best_student.pt',{'model':copy.deepcopy(model.state_dict()),'epoch':epoch,'metrics':metric})
   base.dump(dest/'latest_fidelity.json',{**metric,'gate_passed':passed,'epoch':epoch})
   frame.to_csv(dest/'latest_fidelity_by_slide.csv',index=False)
  history.append(rec);pd.DataFrame(history).to_csv(dest/'training_history.csv',index=False)
  base.save(dest/'resume.pt',{'model':model.state_dict(),'optimizer':opt.state_dict(),'history':history,'epoch':epoch,'best':best,'rng':rng.bit_generator.state})
  print(json.dumps({'tissue':tissue,**rec,'gate_passed':passed}),flush=True)
  if passed:
   base.save(dest/'best_student.pt',{'model':copy.deepcopy(model.state_dict()),'epoch':epoch,'metrics':metric})
   break
 beststate=torch.load(dest/'best_student.pt',map_location=device,weights_only=False);model.load_state_dict(beststate['model'])
 metric,frame=base.evaluate(model,rows,cp,ref,bybatch,batches,device);frame.to_csv(dest/'final_fidelity_by_slide.csv',index=False)
 gate=contract['gate'];passed=metric['rna_jsd']<=gate['rna_jsd_max'] and metric['composition_jsd']<=gate['composition_jsd_max'] and metric['expression_gene_pcc_median']>=gate['expression_gene_pcc_min'] and metric['expression_spot_cosine_median']>=gate['expression_spot_cosine_min'] and metric.get('spatial_theta_gradient_cosine',-1)>=gate['spatial_gradient_cosine_min']
 metric.update(passed=passed,epoch=int(beststate['epoch']),scope=contract['validation_interpretation']);base.dump(dest/'final_fidelity.json',metric)
 # Do not inherit fold-specific scores, stopping epochs or stale teacher provenance.
 keys=['input_dim','hidden','programs','graph_layers','context_graph_layers','graph_heads','type_lift_alpha','sources','technologies','support_mass','hierarchical_logits','support_gamma','neighbors','context_neighbors']
 final={k:cp[k] for k in keys if k in cp};final.update(model={k:v.detach().cpu() for k,v in model.state_dict().items()},inr=None,inr_accepted=False,purpose='same_capacity_ensemble_distilled_deployment_checkpoint',tissue_key=tissue,consensus_reference_sha256=contract['reference_sha256'],gene_order_sha256=contract['gene_order_sha256'],distillation_contract=contract,teacher_fidelity=metric)
 base.save(dest/'model_checkpoint.pt',final);np.save(dest/'program_group_index.npy',groups)
 base.dump(dest/'completed.json',{'tissue':tissue,'publication_gate_passed':passed,'teacher_folds':len(rr),'checkpoint_sha256':base.sha(dest/'model_checkpoint.pt'),'checkpoint_bytes':(dest/'model_checkpoint.pt').stat().st_size,'model_parameter_count':sum(p.numel() for p in model.parameters()),'metrics':metric})
 print(json.dumps({'completed':tissue,'passed':passed,'metrics':metric}),flush=True);CPU.clear();torch.cuda.empty_cache();return True

def main():
 p=argparse.ArgumentParser();p.add_argument('--tissue',action='append');p.add_argument('--epochs',type=int,default=150);p.add_argument('--lr',type=float,default=5e-5);p.add_argument('--seed',type=int,default=20261010);p.add_argument('--output',type=Path,default=OUT/'integrated_ensemble_v1');a=p.parse_args()
 torch.set_num_threads(8);torch.set_float32_matmul_precision('high');torch.cuda.set_per_process_memory_fraction(6*(1<<30)/torch.cuda.get_device_properties(0).total_memory)
 a.output.mkdir(parents=True,exist_ok=True);(a.output/'pid.txt').write_text(str(os.getpid()))
 tissues=a.tissue or sorted(pd.read_csv(OFF/'checkpoints/checkpoint_registry.csv').tissue_key.unique())
 pending=tissues.copy()
 while pending:
  progress=False
  for tissue in pending.copy():
   if run(tissue,a):pending.remove(tissue);progress=True
  if pending and not progress:time.sleep(10)
if __name__=='__main__':main()
