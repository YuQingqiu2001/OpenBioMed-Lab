"""Reference-free coarse RNA to cell-resolution inference adapters.

RNA profiles and continuous axes are fitted to this section's coarse counts.
The generic adapter changes input validation/packing only. BayesGraph fitting,
cross fitting and the verified 1e-9 nullspace projection are frozen upstream.
"""
from pathlib import Path
import json,sys
import h5py,numpy as np,pandas as pd
from scipy import sparse
from .runtime import activate,ENGINE,fresh,sha256

INPUT_SCHEMA='covarst_coarse_counts_and_image_nuclei_v1'
PREPARED_SCHEMA='covarst_coarse_independent_inputs_v1'
CLASSES=['Neoplastic','Inflammatory','Connective','Dead','Epithelial']

def text(h,key):return np.asarray([v.decode() if isinstance(v,bytes) else str(v) for v in h[key][:]])
def read_sparse(h,key,constructor):
    q=h[key]
    value=constructor((q['data'][:].astype(np.float64),q['indices'][:],q['indptr'][:]),shape=tuple(q['shape'][:]))
    value.check_format(full_check=True)
    if np.any(value.data<0) or not np.isfinite(value.data).all():raise ValueError('Nonnegative finite sparse values required')
    return value
def write_sparse(h,key,value):
    q=h.create_group(key);q.attrs['format']=value.format
    for k,v in [('shape',value.shape),('data',value.data),('indices',value.indices),('indptr',value.indptr)]:q.create_dataset(k,data=v,compression='lzf')

def prepare_coarse(input_file,output,sample,seed=20261010):
    activate('cells')
    import prepare_generic_hd_inputs as generic
    core=generic.load_core(ENGINE/'cells')
    with h5py.File(input_file,'r') as h:
        if h.attrs.get('schema')!=INPUT_SCHEMA or int(h.attrs.get('complete',0))!=1:raise ValueError('Expected a complete coarse RNA / real image nucleus input')
        if h.attrs.get('rna_source')!='measured_coarse_st' or bool(h.attrs.get('external_single_cell_reference_used',True)):
            raise ValueError('Only the input ST counts may provide RNA information; external single-cell references are forbidden')
        if not bool(h.attrs.get('real_image_nuclei',False)):raise ValueError('The cell census must come from real image nuclei')
        matrix=read_sparse(h,'matrix',sparse.csc_matrix);geometry=read_sparse(h,'geometry',sparse.csr_matrix)
        genes=text(h,'gene_name');barcodes=text(h,'spots/barcode')
        values={k:h['cells/'+k][:] for k in ['cell_id','class_id','coords_um','spatial_px','features','rna_capacity','owner_bin','domain_component']}
        coords=h['spots/coords_um'][:];components=h['spots/component'][:].astype(np.int32)
    n=len(values['cell_id']);s=len(barcodes);g=len(genes)
    if matrix.shape!=(g,s) or geometry.shape!=(s,n) or n==0 or s<5 or g<20:raise ValueError('Invalid coarse matrix/geometry dimensions')
    if len(np.unique(genes))!=g or len(np.unique(barcodes))!=s or len(np.unique(values['cell_id']))!=n:raise ValueError('Unique aligned gene, spot and cell IDs required')
    if np.any(matrix.data!=np.floor(matrix.data)):raise ValueError('Measured coarse ST requires integer counts, not predicted probabilities')
    library=np.asarray(matrix.sum(0)).ravel()
    if np.any(library<=0):raise ValueError('Provide positive in-tissue coarse spots; record excluded spots upstream')
    for name,shape in [('coords_um',(n,2)),('spatial_px',(n,2)),('rna_capacity',(n,)),('class_id',(n,)),('owner_bin',(n,)),('domain_component',(n,))]:
        if values[name].shape!=shape or not np.isfinite(values[name]).all():raise ValueError('Invalid '+name)
    if coords.shape!=(s,2) or not np.isfinite(coords).all() or components.shape!=(s,) or np.any(components<0):raise ValueError('Invalid physical spot coordinates/components')
    features=values['features']
    if features.ndim!=2 or features.shape[0]!=n or not np.isfinite(features).all():raise ValueError('Invalid standardized real nucleus features')
    classes=values['class_id'].astype(np.int32);base=values['rna_capacity'].astype(np.float64);owner=values['owner_bin'].astype(np.int32)
    if np.any(base<=0) or classes.min()<0 or classes.max()>4 or owner.min()<0 or owner.max()>=s:raise ValueError('Invalid capacity / five-class PanNuke census / owner')
    if not np.array_equal(values['domain_component'],components[owner]):raise ValueError('Cells must be linked within their measured-domain component')
    for spot in range(s):
        cells=geometry.indices[geometry.indptr[spot]:geometry.indptr[spot+1]]
        if np.any(values['domain_component'][cells]!=components[spot]):raise ValueError('Capture geometry crosses measured-domain components')
    _,mass,ncells=core.census_from_geometry(geometry,classes,base);occupied=ncells>0
    if np.count_nonzero(occupied)<5:raise ValueError('Insufficient occupied coarse spots')
    marker_sets={k:sorted(set(v).intersection(genes)) for k,v in generic.DEFAULT_MARKERS.items()}
    total_by_gene=np.asarray(matrix@occupied.astype(np.float64)).ravel();global_mean=total_by_gene/total_by_gene.sum()
    prior=core.marker_prior(global_mean,genes,marker_sets);fit_genes=core.select_fit_genes(matrix,genes,marker_sets)
    # Preserve the existing ~1 mm spatial blocks in physical coordinates.
    positions=pd.DataFrame({'array_col':coords[:,0]/16.,'array_row':coords[:,1]/16.})
    blocks=core.spatial_blocks(positions)
    alpha,profiles,alpha_cv=core.fit_alpha_and_profiles(matrix,library,mass,prior,fit_genes,blocks,.08,5,seed)
    capacity=base*alpha[classes];composition,_,_=core.census_from_geometry(geometry,classes,capacity)
    shrink,shrink_cv=core.hierarchy_shrinkage_cv(matrix,library,composition,prior,fit_genes,blocks,.08,5)
    profiles=core.apply_hierarchy_shrinkage(profiles,shrink,composition,global_mean)
    panel_index=core.select_fit_genes(matrix,genes,marker_sets,limit=min(278,g))
    panel=genes[panel_index]
    summary={'schema':PREPARED_SCHEMA,'sample':sample,'panel_genes':panel.tolist(),'external_reference_used':False,'fine_expression_used':False,'RNA_profile_source':'this coarse ST only','counts_semantics':'measured coarse ST integer counts','coordinate_unit_um':100.,'physical_coordinates_declared':True,'real_image_nuclei':True,'alpha_cv':alpha_cv,'shrinkage_cv':shrink_cv,'unoccupied_spot_mass':float(matrix[:,~occupied].sum()),'class_names':CLASSES,'input_sha256':sha256(input_file),'adapter_sha256':sha256(__file__),'CV_interpretation':'conditional/transductive coarse model selection with global profiles fitted to the same section'}
    output=fresh(output);partial=output.with_name(output.name+'.partial')
    with h5py.File(partial,'x') as h:
        h.attrs.update(schema=PREPARED_SCHEMA,complete=0,coordinate_unit_um=100.,input_provenance_json=json.dumps(summary),summary_json=json.dumps(summary))
        for name,value in [('gene_name',genes),('panel_gene_name',panel),('class_names',CLASSES),('spots/barcode',barcodes)]:h.create_dataset(name,data=np.asarray(value,object),dtype=h5py.string_dtype('utf-8'))
        cells=h.create_group('cells')
        for name,value in {**{k:v for k,v in values.items() if k!='coords_um'},'coords':values['coords_um']/100.,'rna_capacity':capacity,'direct_capture_overlap_fraction':np.asarray(geometry.sum(0)).ravel(),'cross_component_fallback':np.zeros(n,bool)}.items():cells.create_dataset(name,data=value,compression='lzf')
        h.create_dataset('spots/coords',data=coords/100.);h.create_dataset('spots/component',data=components)
        write_sparse(h,'matrix',matrix);write_sparse(h,'geometry',geometry);h.create_dataset('model/class_gene_probability',data=profiles)
        h.attrs['complete']=1
    partial.replace(output)
    return {'output':str(output),'spots':s,'cells':n,'genes':g,'external_reference_used':False}

def inspect_coarse(a):
    """Pack generic inputs into the unchanged original BayesGraph fitter contract."""
    activate('cells');import fit_native16_bayesgraph as native
    with h5py.File(a.prepared,'r') as h:
        if h.attrs.get('schema')!=PREPARED_SCHEMA or int(h.attrs.get('complete',0))!=1:raise ValueError('Expected prepared generic coarse input')
        provenance=json.loads(h.attrs['input_provenance_json'])
        if provenance.get('sample')!=a.sample or provenance.get('external_reference_used') is not False or provenance.get('fine_expression_used') is not False:raise ValueError('Provenance violation')
        a.coordinate_unit_um=100.;a.spatial_buffer=a.spatial_buffer_um/100.
        genes=text(h,'gene_name');panel=text(h,'panel_gene_name');names=text(h,'class_names');lookup={v:i for i,v in enumerate(genes)}
        if names.tolist()!=CLASSES or panel.tolist()!=provenance['panel_genes']:raise ValueError('Class/panel mismatch')
        arrays={name:h[key][:] for name,key in [('cell_id','cells/cell_id'),('class_id','cells/class_id'),('owner_bin','cells/owner_bin'),('coords','cells/coords'),('spot_component','spots/component'),('cell_component','cells/domain_component'),('profiles','model/class_gene_probability')]}
        spotcoords=h['spots/coords'][:];n=len(arrays['cell_id']);s=len(spotcoords)
        if tuple(h['matrix/shape'][:])!=(len(genes),s) or tuple(h['geometry/shape'][:])!=(s,n):raise ValueError('Sparse shape mismatch')
        if not np.allclose(arrays['profiles'].sum(1),1,atol=2e-6):raise ValueError('Coarse-only profiles must sum to one')
        extras={key:h['cells/'+key][:] for key in ['domain_component','direct_capture_overlap_fraction','cross_component_fallback']}
    data={**arrays,'gene_names':genes,'panel_genes':panel,'panel_index':np.asarray([lookup[v] for v in panel],np.int32),'spots':pd.DataFrame({'hex_col':spotcoords[:,0],'hex_row':spotcoords[:,1]}),'cell_extras':extras}
    report={'schema':'covarst_generic_coarse_BayesGraph_adapter_v1','sample':a.sample,'class_names':names.tolist(),'n_cells':n,'n_coarse_spots':s,'n_genes':len(genes),'n_panel_genes':len(panel),'coordinate_unit_um':100.,'preflight_passed':True,'external_reference_used':False,'2um_expression_used_for_fit':False,'input_provenance':provenance,'inputs':{'prepared':native.fingerprint(a.prepared)},'CV_interpretation':provenance['CV_interpretation']}
    return report,data

def fit_coarse(arguments):
    activate('cells')
    import fit_generic_v55_bayesgraph as generic
    import fit_bayesgraph_component_projection as strict
    generic.inspect_inputs=inspect_coarse
    args=list(arguments)
    if args[:1]==['--']:args.pop(0)
    sys.argv=[str(ENGINE/'cells/fit_bayesgraph_component_projection.py'),'--engine','v55']+args
    strict.main()

def read_factors(h):
    if int(h.attrs.get('complete',0))!=1:raise ValueError('Incomplete factor file')
    genes=text(h,'gene_name');classes=h['cells/class_id'][:].astype(np.int32)
    profiles=h['model/class_gene_probability'][:].astype(np.float64)
    loadings=h['model/program_loadings'][:].astype(np.float64)
    count=h['model/program_count_by_class'][:].astype(np.int32)
    if profiles.shape[1]!=len(genes) or loadings.shape[0]!=len(profiles) or loadings.shape[2]!=len(genes) or np.any(profiles<0) or not np.isfinite(profiles).all():raise ValueError('Invalid all-gene factor dimensions')
    return genes,classes,profiles,loadings,count

def factor_logits(state,classes,profiles,loadings,count,start,end):
    logits=np.log(np.maximum(profiles[classes,start:end],1e-12))
    for c in np.unique(classes):
        keep=classes==c;k=int(count[c])
        if k:logits[keep]+=state[keep,:k]@loadings[c,:k,start:end]
    return logits

def materialize(factors,output,cell_block=256,gene_block=512):
    output=fresh(output);partial=output.with_name(output.name+'.partial')
    with h5py.File(factors,'r') as source,h5py.File(partial,'x') as h:
        genes,classes,profiles,loadings,count=read_factors(source);n=len(classes);g=len(genes)
        h.attrs.update(schema='covarst_cell_probability_v1',complete=0,expression_semantics='inferred cell probability; all-gene softmax; not measured counts',factor_sha256=sha256(factors))
        source.copy('cells',h);h.create_dataset('gene_name',data=genes.astype(object),dtype=h5py.string_dtype('utf-8'))
        prob=h.create_dataset('P_post',shape=(n,g),dtype='f4',compression='lzf',chunks=(min(cell_block,n),min(gene_block,g)))
        logs=h.create_dataset('log1p_cp10k',shape=(n,g),dtype='f4',compression='lzf',chunks=prob.chunks)
        max_error=0.
        for i in range(0,n,cell_block):
            state=source['model/cell_program_activity'][i:i+cell_block].astype(np.float64);cls=classes[i:i+cell_block]
            logden=np.full(len(cls),-np.inf)
            for j in range(0,g,gene_block):
                tile=factor_logits(state,cls,profiles,loadings,count,j,j+gene_block)
                from scipy.special import logsumexp
                logden=np.logaddexp(logden,logsumexp(tile,axis=1))
            row_sum=np.zeros(len(cls))
            for j in range(0,g,gene_block):
                tile=np.exp(factor_logits(state,cls,profiles,loadings,count,j,j+gene_block)-logden[:,None]).astype(np.float32)
                prob[i:i+cell_block,j:j+gene_block]=tile;logs[i:i+cell_block,j:j+gene_block]=np.log1p(10000*tile);row_sum+=tile.sum(1,dtype=np.float64)
            max_error=max(max_error,float(np.max(np.abs(row_sum-1))))
        if max_error>2e-6:raise ValueError('All-gene probability normalization failed')
        h.attrs.update(complete=1,max_probability_row_sum_error=max_error)
    partial.replace(output)
    return {'output':str(output),'cells':n,'genes':g,'max_probability_row_sum_error':max_error}

def allocate_mass(prepared,factors,output):
    activate('cells');from cell_resolved_visium55_mv3 import allocate_bin_mv3
    output=fresh(output);partial=output.with_name(output.name+'.partial')
    with h5py.File(prepared,'r') as src,h5py.File(factors,'r') as f,h5py.File(partial,'x') as h:
        if int(src.attrs.get('complete',0))!=1:raise ValueError('Incomplete coarse inputs')
        matrix=read_sparse(src,'matrix',sparse.csc_matrix);geometry=read_sparse(src,'geometry',sparse.csr_matrix)
        genes,classes,profiles,loadings,count=read_factors(f);state=f['model/cell_program_activity'][:]
        if not np.array_equal(src['cells/cell_id'][:],f['cells/cell_id'][:]) or not np.array_equal(text(src,'gene_name'),genes):raise ValueError('Prepared counts and factors must share exact cell/gene order')
        capacity=src['cells/rna_capacity'][:];n=len(classes);g=len(genes)
        h.attrs.update(schema='covarst_conserved_cell_mass_v1',complete=0,expression_semantics='fractional captured counts; conservation accounting only')
        f.copy('cells',h);h.create_dataset('gene_name',data=genes.astype(object),dtype=h5py.string_dtype('utf-8'))
        values=h.create_dataset('X_mass',shape=(n,g),dtype='f8',chunks=(min(128,n),min(512,g)),compression='lzf')
        # Preserve the per-observation ledger: a cell may overlap multiple spots.
        links=h.create_group('spot_cell_gene_allocation');unassigned=np.zeros(g);max_error=0.
        for spot in range(matrix.shape[1]):
            cells,gidx,counts,q=allocate_bin_mv3(matrix,spot,geometry,capacity,classes,state,profiles,loadings,count,.25,4)
            if not len(gidx):continue
            if not len(cells):unassigned[gidx]+=counts;continue
            mass=q.astype(np.float64)*counts[None,:]
            error=np.max(np.abs(mass.sum(0)-counts)/np.maximum(counts,1));max_error=max(max_error,float(error))
            if error>2e-6:raise ValueError('Spot-gene conservation failed')
            link=links.create_group(str(spot))
            for key,v in [('cell_index',cells),('gene_index',gidx),('allocated_counts',mass)]:link.create_dataset(key,data=v,compression='lzf')
            for local,cell in enumerate(cells):
                old=values[int(cell),:];old[gidx]+=mass[local];values[int(cell),:]=old
        h.create_dataset('unassigned_gene_mass',data=unassigned);h.attrs.update(complete=1,max_spot_gene_relative_error=max_error,unassigned_total_mass=float(unassigned.sum()))
    partial.replace(output)
    return {'output':str(output),'max_spot_gene_relative_error':max_error,'unassigned_total_mass':float(unassigned.sum())}
