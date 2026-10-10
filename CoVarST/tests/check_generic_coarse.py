"""Small full generic adapter smoke test, explicitly synthetic software data."""
from pathlib import Path
import sys,tempfile,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import covarst.cells as cells
import numpy as np,h5py
from scipy import sparse
import torch
torch.set_num_threads(8)
with tempfile.TemporaryDirectory(prefix='covarst_generic_smoke_') as tmp:
 tmp=Path(tmp);raw=tmp/'raw.h5';prepared=tmp/'prepared.h5';out=tmp/'fit';rng=np.random.default_rng(13);ns=45;nc=ns*5;ng=30
 positions=np.array([(i%9*100,i//9*100) for i in range(ns)],float)
 matrix=sparse.csc_matrix(rng.poisson(8,size=(ng,ns)));geometry=sparse.csr_matrix((np.ones(nc),(np.repeat(np.arange(ns),5),np.arange(nc))),shape=(ns,nc))
 with h5py.File(raw,'w') as h:
  h.attrs.update(schema=cells.INPUT_SCHEMA,complete=1,rna_source='measured_coarse_st',external_single_cell_reference_used=False,real_image_nuclei=True)
  cells.write_sparse(h,'matrix',matrix);cells.write_sparse(h,'geometry',geometry)
  h.create_dataset('gene_name',data=np.asarray([f'SOFTWARE_G{i}' for i in range(ng)],object),dtype=h5py.string_dtype())
  h.create_dataset('spots/barcode',data=np.asarray([f'SOFTWARE_S{i}' for i in range(ns)],object),dtype=h5py.string_dtype());h.create_dataset('spots/coords_um',data=positions);h.create_dataset('spots/component',data=np.zeros(ns,int))
  for name,value in [('cell_id',np.arange(nc)),('class_id',np.tile(np.arange(5),ns)),('coords_um',np.repeat(positions,5,axis=0)+rng.normal(size=(nc,2))),('spatial_px',np.repeat(positions,5,axis=0)),('features',rng.normal(size=(nc,4))),('rna_capacity',np.ones(nc)),('owner_bin',np.repeat(np.arange(ns),5)),('domain_component',np.zeros(nc,int))]:h.create_dataset('cells/'+name,data=value)
 cells.prepare_coarse(raw,prepared,'SOFTWARE_SMOKE',13)
 cells.fit_coarse(['--sample','SOFTWARE_SMOKE','--prepared',str(prepared),'--output-dir',str(out),'--device','cpu','--ranks','0,8','--trend-lambdas','0.05','--spatial-folds','2','--cv-steps-v55','3','--final-steps-v55','3','--fit-batch-size','8','--trend-batch-size','16','--permutations','2','--spatial-buffer-um','0'])
 factors=list(out.glob('*.factorized.h5'));assert len(factors)==2
 for index,factor in enumerate(factors):
  cells.materialize(factor,tmp/f'post{index}.h5');cells.allocate_mass(prepared,factor,tmp/f'mass{index}.h5')
 (ROOT/'validation/generic_coarse_smoke.json').write_text(json.dumps({'status':'passed','data':'synthetic software fixture only; no biological validation claim','workflow':'prepare-coarse -> frozen fitter ranks0,8 -> P_post -> X_mass','short_fit_steps':3,'validation_scope':'entry point, dimension, numerical and output-semantic integration; not learned accuracy'},indent=2)+'\n')
 print('PASS: full generic-coarse packaging workflow, rank0 and rank8 software fixtures')
