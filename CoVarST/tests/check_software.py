"""Meaningful release checks; synthetic fixtures are software tests only."""
from pathlib import Path
import json,sys,tempfile,compileall,importlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from covarst.runtime import activate
import covarst.cells as cells
import numpy as np,h5py
from scipy import sparse

def main():
 assert compileall.compile_dir(str(ROOT/'src'),quiet=1)
 activate('he','reference','reference/upstream_pipeline','cells')
 for name in ['train_virchow2_pan_atlas_consensus_mapper','run_virchow2_bayestme_global_mapper','run_virchow2_fseg_k64_mapper','fit_generic_v55_bayesgraph','fit_native16_bayesgraph','component_nullspace_projector','virchow2_feature_utils','wsi_grid_math']:importlib.import_module(name)
 with tempfile.TemporaryDirectory(prefix='covarst_software_test_') as td:
  td=Path(td);prepared=td/'coarse.h5';factor=td/'factor.h5';massout=td/'mass.h5';post=td/'post.h5'
  counts=np.array([[10,0,7],[3,11,2],[1,2,5]],float)
  mat=sparse.csc_matrix(counts);geom=sparse.csr_matrix(np.array([[1,1,0],[0,0,1],[0,0,0]],float))
  profiles=np.array([[.7,.2,.1],[.1,.7,.2]],np.float32);cls=np.array([0,1,1]);ids=np.arange(3);genes=np.array(['TEST_A','TEST_B','TEST_C'])
  with h5py.File(prepared,'w') as h:
   h.attrs['complete']=1;cells.write_sparse(h,'matrix',mat);cells.write_sparse(h,'geometry',geom)
   h.create_dataset('cells/cell_id',data=ids);h.create_dataset('cells/rna_capacity',data=np.array([1,2,1]));h.create_dataset('gene_name',data=genes.astype(object),dtype=h5py.string_dtype())
  with h5py.File(factor,'w') as h:
   h.attrs['complete']=1;h.create_dataset('gene_name',data=genes.astype(object),dtype=h5py.string_dtype())
   for k,v in [('cell_id',ids),('class_id',cls),('spatial',np.zeros((3,2)))]:h.create_dataset('cells/'+k,data=v)
   for k,v in [('class_gene_probability',profiles),('program_loadings',np.zeros((2,1,3))),('program_count_by_class',np.array([0,0])),('cell_program_activity',np.zeros((3,1)))]:h.create_dataset('model/'+k,data=v)
  result=cells.materialize(factor,post,2,2)
  with h5py.File(post) as h:assert np.allclose(h['P_post'][:],profiles[cls]);assert result['max_probability_row_sum_error']<1e-6
  result=cells.allocate_mass(prepared,factor,massout)
  with h5py.File(massout) as h:
   assert np.allclose(h['X_mass'][:].sum(0)+h['unassigned_gene_mass'][:],counts.sum(1),atol=1e-6)
   assert np.allclose(h['unassigned_gene_mass'][:],counts[:,2])
  import component_nullspace_projector as projector
  projector.self_test()
  # Generic-coarse preparation and metadata preflight: no external RNA field.
  rng=np.random.default_rng(7);raw=td/'generic_coarse.h5';ready=td/'prepared.h5';ng=30;ns=25;nc=125
  matrix=sparse.csc_matrix(rng.poisson(6,size=(ng,ns)))
  geometry=sparse.csr_matrix((np.ones(nc),(np.repeat(np.arange(ns),5),np.arange(nc))),shape=(ns,nc))
  spots=np.array([(i%5*100,i//5*100) for i in range(ns)],float)
  with h5py.File(raw,'w') as h:
   h.attrs.update(schema=cells.INPUT_SCHEMA,complete=1,rna_source='measured_coarse_st',external_single_cell_reference_used=False,real_image_nuclei=True)
   cells.write_sparse(h,'matrix',matrix);cells.write_sparse(h,'geometry',geometry)
   h.create_dataset('gene_name',data=np.asarray([f'SOFTWARE_TEST_G{i}' for i in range(ng)],object),dtype=h5py.string_dtype())
   h.create_dataset('spots/barcode',data=np.asarray([f'SOFTWARE_TEST_S{i}' for i in range(ns)],object),dtype=h5py.string_dtype());h.create_dataset('spots/coords_um',data=spots);h.create_dataset('spots/component',data=np.zeros(ns,int))
   for k,v in [('cell_id',np.arange(nc)),('class_id',np.tile(np.arange(5),ns)),('coords_um',np.repeat(spots,5,axis=0)+rng.normal(size=(nc,2))),('spatial_px',np.repeat(spots,5,axis=0)),('features',rng.normal(size=(nc,4))),('rna_capacity',np.ones(nc)),('owner_bin',np.repeat(np.arange(ns),5)),('domain_component',np.zeros(nc,int))]:h.create_dataset('cells/'+k,data=v)
  cells.prepare_coarse(raw,ready,'SOFTWARE_TEST',seed=7)
  import argparse
  report,data=cells.inspect_coarse(argparse.Namespace(prepared=ready,sample='SOFTWARE_TEST',spatial_buffer_um=150.))
  assert report['n_cells']==nc and report['external_reference_used'] is False
  with h5py.File(raw,'r+') as h:h.attrs['external_single_cell_reference_used']=True
  try:cells.prepare_coarse(raw,td/'forbidden.h5','SOFTWARE_TEST')
  except ValueError:pass
  else:raise AssertionError('External RNA reference was not rejected')
  import inspect
  signatures={'projector':str(inspect.signature(projector.measurement_nullspace_center))}
 report={'status':'passed','test_data':'explicit synthetic software fixtures only, not study figures or biological validation','source_compile':True,'frozen_engine_imports':True,'P_post_all_gene_normalization':True,'X_mass_spot_gene_conservation_and_unassigned_mass':True,'strict_nullspace_oracle':True,'generic_coarse_preparation_contract':True,'external_RNA_reference_rejected':True,**signatures}
 (ROOT/'validation/software_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
