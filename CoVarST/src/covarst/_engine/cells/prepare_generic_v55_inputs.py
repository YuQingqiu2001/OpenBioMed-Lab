#!/usr/bin/env python3
"""Prepare independent virtual55 coarse inputs from native16 and real nuclei.

No fine expression input exists. Original circle integration, Voronoi geometry,
and global-profile functions are reused without editing upstream source.
Cells between 55um circles are retained and explicitly flagged as interpolation.
Run only in existing ST_GJH. New output directory; no overwrite.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='8'
if hasattr(os,'sched_setaffinity'):
    os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
import argparse
import importlib
import json
from pathlib import Path
import sys
import time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import prepare_generic_hd_inputs as generic

SCHEMA='generic_hd_virtual55_independent_inputs_v1'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',required=True)
    for key in ('matrix','positions','scalefactors','features','output-dir'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--pipeline-dir',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--coverage-threshold',type=float,default=.95)
    p.add_argument('--embedding-pcs',type=int,default=16)
    p.add_argument('--pca-sample-per-class',type=int,default=5000)
    p.add_argument('--class-prior-strength',type=float,default=.08)
    p.add_argument('--cv-folds',type=int,default=5)
    p.add_argument('--panel-size',type=int,default=278)
    p.add_argument('--seed',type=int,default=20260828)
    a=p.parse_args();a.mode='coarse'
    if not 0<a.coverage_threshold<=1:
        p.error('Invalid circle coverage gate')
    for key in ('matrix','positions','scalefactors','features'):
        source=getattr(a,key)
        if not source.is_file():raise FileNotFoundError(source)
        if any(token in str(source).lower() for token in ('002um','direct2um','calibrat')):
            raise ValueError(f'Fine-expression/calibrated source forbidden: {source}')
    if a.output_dir.exists():raise FileExistsError(a.output_dir)
    a.output_dir.mkdir(parents=True)
    started=time.time()
    factors,_=generic.check_input_resolution(a,16)
    mpp=float(factors['microns_per_pixel'])
    all_ids,xy,raw_classes=generic.feature_metadata(a.features,mpp)
    core=generic.load_core(a.pipeline_dir)
    sys.path.insert(0,str(a.pipeline_dir))
    virtual=importlib.import_module('cell_resolved_visium55')
    virtual.core=core
    matrix,library,barcodes,gene_ids,genes,types,source_columns,positions,domain_audit=generic.load_positive_tissue(a,core)
    _,_,physical_affine=generic.affine_from_positions(positions,mpp,16)
    position_path=a.output_dir/f'{a.sample}.source16.positions.tsv.gz'
    positions.to_csv(position_path,sep='\t',index=False,compression='gzip')
    positions,owner16,eligible,uv_all,affine,spatial_audit=core.align_positions(position_path,barcodes,xy)
    eligible &= raw_classes>0
    if not eligible.any():raise ValueError('No eligible real nuclei')
    print(f'[{a.sample}] Creating virtual55 circles from {matrix.shape[1]:,} raw16 bins',flush=True)
    counts,source_overlap,spots,gap_gene,degradation=virtual.build_virtual_spots(matrix,positions,affine,a.sample,a.coverage_threshold)
    source_component=core.measured_components(positions)
    spot_component=virtual.spot_source_components(source_overlap,source_component)
    cell_component_all=np.full(len(all_ids),-1,np.int32)
    cell_component_all[eligible]=source_component[owner16[eligible]]
    centers_uv=spots[['array_col','array_row']].to_numpy(np.float64)
    owner55,cross=virtual.nearest_spot_within_component(uv_all,eligible,cell_component_all,centers_uv,spot_component)
    nuclei=generic.load_nuclei(a.features,owner55,eligible,uv_all)
    nuclei['domain_component']=cell_component_all[eligible]
    print(f'[{a.sample}] Constructing measured-domain Voronoi geometry for {eligible.sum():,} nuclei',flush=True)
    geometry,geometry_audit=virtual.build_full_domain_voronoi_overlap(nuclei,affine,centers_uv,spot_component,cross,55/2/16)
    base=core.base_capacity(nuclei)
    _,mass,n_cells=core.census_from_geometry(geometry,nuclei['class_id'],base)
    occupied=n_cells>0
    coarse_library=np.asarray(counts.sum(axis=0)).ravel()
    marker_sets={name:sorted(set(values).intersection(genes)) for name,values in generic.DEFAULT_MARKERS.items()}
    total_by_gene=np.asarray(counts@occupied.astype(np.float64)).ravel()
    global_mean=total_by_gene/total_by_gene.sum()
    prior=core.marker_prior(global_mean,genes,marker_sets)
    fit_genes=core.select_fit_genes(counts,genes,marker_sets)
    blocks=core.spatial_blocks(spots)
    print(f'[{a.sample}] Independently fitting global profiles from {counts.shape[1]:,} virtual55 spots',flush=True)
    alpha,profiles,alpha_cv=core.fit_alpha_and_profiles(counts,coarse_library,mass,prior,fit_genes,blocks,a.class_prior_strength,a.cv_folds,a.seed)
    capacity=(base*alpha[nuclei['class_id']]).astype(np.float32)
    composition,_,_=core.census_from_geometry(geometry,nuclei['class_id'],capacity)
    shrinkage,shrinkage_cv=core.hierarchy_shrinkage_cv(counts,coarse_library,composition,prior,fit_genes,blocks,a.class_prior_strength,a.cv_folds)
    profiles=core.apply_hierarchy_shrinkage(profiles,shrinkage,composition,global_mean)
    if not np.isfinite(profiles).all() or np.any(profiles<0) or not np.allclose(profiles.sum(1),1,atol=2e-6):
        raise ValueError('Invalid independently fitted class probabilities')
    features,_,_=core.learn_cell_features(a.features,nuclei,a.embedding_pcs,a.pca_sample_per_class,a.seed)
    if not np.isfinite(features).all():raise ValueError('Invalid image features')
    indices=core.select_fit_genes(counts,genes,marker_sets,limit=a.panel_size)
    panel=list(dict.fromkeys(genes[indices].tolist()))
    if len(panel)<20:raise ValueError('Insufficient unique coarse-only panel genes')
    input_provenance={name:generic.fingerprint(getattr(a,name)) for name in ('matrix','positions','scalefactors','features')}
    source_provenance={name:generic.fingerprint(a.pipeline_dir/name) for name in ('cell_resolved_v3.py','cell_resolved_visium55.py')}
    summary={'schema':SCHEMA,'sample':a.sample,'complete':True,'class_names':generic.CLASSES,
        'source_resolution_um':16,'virtual_spot_diameter_um':55,'virtual_spot_pitch_um':100,
        'counts_semantics':'fractional captured UMI mass: exact circle/16um-square area integration; not native integer55 measurement',
        'profile_source':'independent fit of virtual55 measurements; no fitted native16 profiles/programs/expression reused',
        'source_domain':domain_audit,'spatial_audit':spatial_audit,'physical_affine':physical_affine,
        'degradation':degradation,'geometry':geometry_audit,'n_cells':len(nuclei['cell_id']),
        'n_genes':len(genes),'n_virtual_spots':len(spots),'n_unoccupied_spots':int((~occupied).sum()),
        'unassigned_mass_unoccupied_spots':float(counts[:,~occupied].sum()),'marker_sets':marker_sets,
        'alpha':alpha.tolist(),'alpha_cv':alpha_cv,'shrinkage':shrinkage.tolist(),'shrinkage_cv':shrinkage_cv,
        'panel_genes':panel,'panel_selection':'virtual55 coarse abundance plus explicit broad lineage anchors; no fine expression',
        'coordinate_unit_um':100,'coordinate_semantics':'physical Euclidean array_col/row times16/100; not hex lattice indices',
        'microns_per_pixel':mpp,'input_provenance':input_provenance,'original_sources':source_provenance,
        'adapter_sha256':generic.digest(__file__),'2um_expression_used_for_fit':False,
        'CV_interpretation':'conditional within-section coarse-profile/model selection; not independent single-cell validation',
        'elapsed_seconds':time.time()-started,'seed':a.seed}
    target=a.output_dir/f'{a.sample}.virtual55.prepared.h5'
    partial=target.with_name(target.name+'.partial')
    text=h5py.string_dtype('utf-8')
    with h5py.File(partial,'x') as h:
        h.attrs.update(schema=SCHEMA,complete=0,coordinate_unit_um=100.,microns_per_pixel=mpp,
            input_provenance_json=json.dumps(summary,ensure_ascii=False),summary_json=json.dumps(summary,ensure_ascii=False))
        h.create_dataset('gene_name',data=np.asarray(genes,object),dtype=text)
        h.create_dataset('class_names',data=np.asarray(generic.CLASSES,object),dtype=text)
        h.create_dataset('panel_gene_name',data=np.asarray(panel,object),dtype=text)
        cells=h.create_group('cells')
        for name,value in {'cell_id':nuclei['cell_id'],'source_cell_index':nuclei['source_index'],
            'class_id':nuclei['class_id'],'owner_bin':nuclei['owner_bin'],'rna_capacity':capacity,
            'spatial_px':nuclei['centroid_xy'],'coords':nuclei['centroid_uv']*16/100,
            'features':features,'native16_array_uv':nuclei['centroid_uv'],
            'domain_component':nuclei['domain_component'],
            'direct_capture_overlap_fraction':nuclei['direct_capture_overlap_fraction_sum'],
            'cross_component_fallback':nuclei['spatial_evidence_level']==2,
            'spatial_evidence_level':nuclei['spatial_evidence_level']}.items():
            cells.create_dataset(name,data=value,compression='lzf')
        group=h.create_group('spots');group.create_dataset('coords',data=centers_uv*16/100)
        group.create_dataset('native16_array_colrow',data=centers_uv)
        group.create_dataset('component',data=spot_component)
        group.create_dataset('barcode',data=np.asarray(spots.barcode,object),dtype=text)
        for name,mat in [('matrix',counts),('geometry',geometry),('source16_overlap',source_overlap)]:
            group=h.create_group(name)
            for key,value in {'shape':mat.shape,'data':mat.data,'indices':mat.indices,'indptr':mat.indptr}.items():
                group.create_dataset(key,data=value,compression='lzf')
            group.attrs['format']='csc' if name=='matrix' else 'csr'
        h.create_dataset('model/class_gene_probability',data=profiles)
        h.create_dataset('degradation/source_positive_raw16_column_index',data=source_columns,compression='lzf')
        h.create_dataset('degradation/gap_gene_mass',data=gap_gene)
        h.attrs['complete']=1
    partial.replace(target)
    spots.to_csv(a.output_dir/f'{a.sample}.virtual55.positions.tsv.gz',sep='\t',index=False,compression='gzip')
    summary['output']=str(target)
    generic.write_json(a.output_dir/f'{a.sample}.virtual55.preparation.qc.json',summary)
    print(json.dumps({'complete':True,'output':str(target),'n_cells':len(nuclei['cell_id']),'n_spots':len(spots),'elapsed_seconds':time.time()-started}),flush=True)


if __name__=='__main__':main()
