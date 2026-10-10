"""Align frozen cross-fitted teacher outputs with unchanged H&E features."""
from pathlib import Path
import csv, hashlib, json, time, sys, os
import numpy as np, pandas as pd
import os
from pathlib import Path
ROOT=Path(os.environ['COVARST_STUDY_ROOT']).resolve()
OFF=ROOT
OUT=Path(os.environ['COVARST_BUILD_ROOT']).resolve()



from run_virchow2_bayestme_global_mapper import _morphology_graph
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 rows=pd.read_csv(OFF/'checkpoints/checkpoint_registry.csv')
 OUT.mkdir(parents=True,exist_ok=True); manifest=[]
 for tissue,group in rows.groupby('tissue_key',sort=True):
  dest=OUT/'teacher_cache'/tissue;dest.mkdir(parents=True,exist_ok=True)
  seen=set(); records=[]
  for row in group.itertuples(index=False):
   fold=OFF/str(row.checkpoint_relative_path).replace('\\','/');fold=fold.parent.parent
   complete=json.loads((fold/'fold_complete.json').read_text())
   assert complete['status']=='complete'
   assert complete['fixed_reference_sha256']==row.fixed_reference_sha256
   inference=fold/'pure_inference_test'
   info=json.loads((inference/'inference_manifest.json').read_text())
   assert info['transcriptomic_inputs_loaded'] is False
   assert info['checkpoint_sha256']==row.checkpoint_sha256
   table=pd.read_csv(inference/'inference_spot_index.csv')
   slides=pd.read_csv(inference/'inference_slides.csv')
   theta=np.load(inference/'theta_rna_predicted.npy',mmap_mode='r')
   composition=np.load(inference/'theta_composition_predicted.npy',mmap_mode='r')
   assert theta.shape==composition.shape==(len(table),int(row.programs))
   assert set(table.patient.astype(str))=={str(row.heldout_patient)}
   for sr in slides.itertuples(index=False):
    key=str(sr.slide_key);assert key not in seen;seen.add(key)
    block=table.loc[table.inference_slide_key==key].sort_values('prediction_row')
    tag=hashlib.sha256(key.encode()).hexdigest()[:16]
    cache=dest/f'{tag}.npz'; meta=cache.with_suffix('.json')
    provenance={'tissue':tissue,'slide':key,'patient':str(sr.patient),'source':str(sr.source),'technology':str(block.technology.iloc[0]),'spots':len(block),'programs':int(row.programs),'teacher_fold':int(row.fold_index),'teacher_checkpoint_sha256':str(row.checkpoint_sha256),'reference_sha256':str(row.fixed_reference_sha256),'teacher_patient_held_out':True,'teacher_inr_accepted':str(row.inr_accepted).lower()=='true','feature_path':str(sr.feature_path),'teacher_rna_sha256':sha(inference/'theta_rna_predicted.npy'),'teacher_composition_sha256':sha(inference/'theta_composition_predicted.npy'),'index_sha256':sha(inference/'inference_spot_index.csv')}
    if cache.exists() and meta.exists():
     old=json.loads(meta.read_text())
     assert all(old[k]==v for k,v in provenance.items()),f'stale teacher cache: {key}'
     records.append(old);continue
    path=Path(sr.feature_path)
    with np.load(path,allow_pickle=False) as z:
     barcodes=z['barcode'].astype(str);lookup={v:i for i,v in enumerate(barcodes)}
     ix=np.array([lookup[v] for v in block.barcode.astype(str)],np.int64)
     f=np.array(z['features'][ix],np.float16);coords=np.array(z['coords'][ix],np.float32)
    r=np.array(theta[block.prediction_row],np.float32);c=np.array(composition[block.prediction_row],np.float32)
    assert f.shape==(len(block),2560) and np.isfinite(f).all()
    assert np.isfinite(r).all() and np.isfinite(c).all() and (r>=0).all() and (c>=0).all()
    assert np.allclose(r.sum(1),1,atol=2e-6) and np.allclose(c.sum(1),1,atol=2e-6)
    local=_morphology_graph(coords,f.astype(np.float32),neighbors=6)
    context=_morphology_graph(coords,f.astype(np.float32),neighbors=18)
    arrays=dict(features=f,coords=coords,teacher_rna=r,teacher_composition=c,barcode=block.barcode.astype(str).to_numpy(dtype='U'),global_index=block.global_index.to_numpy(np.int64))
    for prefix,graph in [('local',local),('context',context)]:
     for label,value in zip(['source','target','prior'],graph):arrays[f'{prefix}_{label}']=value
    with cache.with_suffix('.partial').open('wb') as h:np.savez(h,**arrays)
    cache.with_suffix('.partial').replace(cache)
    provenance.update(cache=str(cache),cache_sha256=sha(cache),feature_sha256=sha(path))
    meta.write_text(json.dumps(provenance,indent=2,ensure_ascii=False))
    records.append(provenance)
    print(json.dumps({'cached':key,'tissue':tissue,'spots':len(block),'teacher_fold':int(row.fold_index)}),flush=True)
   del theta,composition
  manifest.extend(records)
  (dest/'manifest.json').write_text(json.dumps(records,indent=2,ensure_ascii=False))
  print(json.dumps({'tissue_complete':tissue,'slides':len(records),'patients':len({r['patient'] for r in records}),'spots':sum(r['spots'] for r in records),'largest_slide':max(r['spots'] for r in records)}),flush=True)
 (OUT/'teacher_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
 pd.DataFrame([{k:v for k,v in r.items() if k not in ('feature_path','cache')} for r in manifest]).to_csv(OUT/'teacher_alignment_audit.csv',index=False)
 print(json.dumps({'prepared_slides':len(manifest),'spots':sum(r['spots'] for r in manifest),'teacher_folds':len(rows)}),flush=True)
if __name__=='__main__':main()
