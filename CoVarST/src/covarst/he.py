"""Cancer-specific H&E-only inference and fixed-reference expression decoding."""
from pathlib import Path
import json
import numpy as np
from .runtime import activate,asset_root,sha256,fresh

def load_model(cancer,assets=None,device='cpu'):
    import torch
    activate('he')
    from train_virchow2_pan_atlas_consensus_mapper import PanAtlasDualMapper
    root=asset_root(assets);folder=root/'models'/cancer
    cp=torch.load(folder/'model_checkpoint.pt',map_location='cpu',weights_only=False)
    if cp.get('purpose')!='same_capacity_ensemble_distilled_deployment_checkpoint' or cp.get('tissue_key')!=cancer:
        raise ValueError('Checkpoint cancer or deployment purpose mismatch')
    if not cp.get('teacher_fidelity',{}).get('passed',False):raise ValueError('Checkpoint has not passed its declared consolidation fidelity gate')
    refroot=root/'references'/cancer
    for name,key in [('consensus_reference_probability_10k.npy','consensus_reference_sha256'),('gene_order_10k.npy','gene_order_sha256')]:
        if sha256(refroot/name)!=cp[key]:raise ValueError(f'Frozen reference mismatch: {name}')
    groups=np.load(folder/'program_group_index.npy',allow_pickle=False)
    model=PanAtlasDualMapper(int(cp['input_dim']),hidden=int(cp['hidden']),programs=int(cp['programs']),program_groups=groups,graph_layers=int(cp['graph_layers']),context_graph_layers=int(cp['context_graph_layers']),graph_heads=int(cp['graph_heads']),dropout=.10,type_lift_alpha=float(cp['type_lift_alpha']),sources=int(cp['sources']),technologies=int(cp['technologies']),support_mass=float(cp['support_mass']),hierarchical_logits=bool(cp.get('hierarchical_logits',False)))
    model.load_state_dict(cp['model'],strict=True);model.eval().to(device)
    return model,cp,refroot

def predict_programs(features,coords,model,cp,device='cpu'):
    import torch
    activate('he')
    from run_virchow2_bayestme_global_mapper import _morphology_graph
    from train_virchow2_pan_atlas_consensus_mapper import _to_device_graph_pair
    x=np.asarray(features,np.float32);xy=np.asarray(coords,np.float32)
    if x.ndim!=2 or x.shape[1]!=int(cp['input_dim']) or xy.shape!=(len(x),2) or len(x)<2:
        raise ValueError('Expected at least two spots, features N x 2560, coordinates N x 2')
    if not np.isfinite(x).all() or not np.isfinite(xy).all():raise ValueError('Nonfinite H&E inputs')
    graphs=tuple(_morphology_graph(xy,x,neighbors=int(cp[key])) for key in ['neighbors','context_neighbors'])
    with torch.inference_mode():
        r,c,_,aux=model(torch.as_tensor(x,device=device),_to_device_graph_pair(graphs,device),support_gamma=float(cp.get('support_gamma',1.)),domain_strength=0.)
        result=[torch.softmax(r,1).cpu().numpy(),torch.softmax(c,1).cpu().numpy(),torch.sigmoid(aux['support_logits']).cpu().numpy()]
    if not all(np.isfinite(v).all() for v in result):raise ValueError('Nonfinite model output')
    return result

def infer_spots(cancer,feature_file,output,assets=None,device='cpu',batch=None,gene_block=512,spot_block=2048):
    import h5py,pandas as pd
    model,cp,refroot=load_model(cancer,assets,device)
    with np.load(feature_file,allow_pickle=False) as z:
        f=z['features'];xy=z['coords'];barcodes=z['barcode'].astype(str)
    if len(barcodes)!=len(f) or len(np.unique(barcodes))!=len(f):raise ValueError('Unique barcodes must align to features')
    r,c,s=predict_programs(f,xy,model,cp,device)
    w=np.load(refroot/'consensus_reference_probability_10k.npy',mmap_mode='r')
    if batch is not None:
        table=pd.read_csv(refroot/'batch_adapter_metadata.csv');rows=table[table.batch.astype(str)==batch]
        if len(rows)!=1:raise ValueError('Unknown source|technology adapter; omit --batch for the neutral reference')
        w=np.load(refroot/'consensus_reference_probability_by_batch_10k.npy',mmap_mode='r')[int(rows.iloc[0].batch_index)]
    genes=np.load(refroot/'gene_order_10k.npy',allow_pickle=False).astype(str)
    if w.shape!=(r.shape[1],len(genes)) or not np.allclose(w.sum(1),1,atol=3e-6):raise ValueError('Invalid fixed reference probabilities')
    output=fresh(output);partial=output.with_name(output.name+'.partial')
    with h5py.File(partial,'x') as h:
        h.attrs.update(schema='covarst_spot_expression_v1',complete=0,cancer=cancer,expression_semantics='fixed-reference expression probabilities and log1p CP10K; not measured counts',transcriptomic_inputs_loaded=False,reference_sha256=cp['consensus_reference_sha256'],gene_order_sha256=cp['gene_order_sha256'],adapter=batch or 'neutral')
        for name,v in [('theta_rna',r),('theta_composition',c),('program_support_probability',s),('spatial',xy)]:h.create_dataset(name,data=v,compression='lzf')
        for name,v in [('barcode',barcodes),('gene_name',genes)]:h.create_dataset(name,data=v.astype(object),dtype=h5py.string_dtype('utf-8'))
        chunks=(min(spot_block,len(f)),min(gene_block,len(genes)))
        prob=h.create_dataset('expression_probability',shape=(len(f),len(genes)),dtype='f4',chunks=chunks,compression='lzf');log=h.create_dataset('expression_log1p_cp10k',shape=prob.shape,dtype='f4',chunks=chunks,compression='lzf')
        for start in range(0,len(f),spot_block):
            for j in range(0,len(genes),gene_block):
                tile=r[start:start+spot_block]@np.asarray(w[:,j:j+gene_block],np.float32)
                prob[start:start+spot_block,j:j+gene_block]=tile;log[start:start+spot_block,j:j+gene_block]=np.log1p(10000*tile)
        h.attrs['complete']=1
    partial.replace(output)
    return {'output':str(output),'cancer':cancer,'spots':len(f),'genes':len(genes),'programs':r.shape[1],'transcriptomic_inputs_loaded':False}

def extract_wsi(wsi,output,paramnet_root,virchow2_root,device='cuda',mpp=None,batch_size=16):
    import torch,openslide
    from PIL import Image
    activate('he')
    from wsi_grid_math import segment_tissue,build_patch_grid,slide_mpp
    from extract_hest_breast_paramnet_uni_features import load_paramnet,normalize_batch
    from virchow2_feature_utils import load_virchow2_model,virchow2_features
    output=fresh(output);dev=torch.device(device)
    paramnet=load_paramnet(Path(paramnet_root),dev);encoder,meta=load_virchow2_model(virchow2_root,dev)
    mean=torch.tensor([.485,.456,.406],device=dev)[None,:,None,None];std=torch.tensor([.229,.224,.225],device=dev)[None,:,None,None]
    with openslide.OpenSlide(str(wsi)) as slide:
        width,height=slide.dimensions
        if mpp is None:mx,my,mpp_source=slide_mpp(slide)
        else:
            if not np.isfinite(mpp) or mpp<=0:raise ValueError('Positive measured MPP required')
            mx=my=float(mpp);mpp_source='explicit_measured_mpp'
        thumb=np.asarray(slide.get_thumbnail((2048,2048)).convert('RGB'));mask,mask_audit=segment_tissue(thumb)
        grid,grid_audit=build_patch_grid(width,height,mx,my,mask,thumb.shape[1]/width,thumb.shape[0]/height,112.,100.,.20)
        features=[]
        for start in range(0,len(grid),batch_size):
            images=[]
            for row in grid.iloc[start:start+batch_size].itertuples():
                patch=slide.read_region((int(row.patch_left_x),int(row.patch_left_y)),0,(int(row.patch_width_src),int(row.patch_height_src))).convert('RGB').resize((224,224),Image.Resampling.LANCZOS)
                images.append(np.asarray(patch))
            images=normalize_batch(paramnet,np.stack(images),dev)
            features.append(virchow2_features(encoder,images,dev,mean,std))
            print(f'WSI features: {min(start+batch_size,len(grid))}/{len(grid)}',flush=True)
    coords=grid[['patch_center_x','patch_center_y']].to_numpy(np.float32);ids=np.asarray([f'virtual_spot_{i:07d}' for i in range(len(grid))])
    with output.open('xb') as stream:np.savez_compressed(stream,features=np.concatenate(features),coords=coords,barcode=ids,coords_um=coords*np.array([mx,my]))
    grid.to_csv(output.with_suffix('.grid.csv'),index=False)
    output.with_suffix('.metadata.json').write_text(json.dumps({'mpp':[mx,my],'mpp_source':mpp_source,'feature_model':meta,'grid':grid_audit,'tissue':mask_audit,'RNA_inputs':False},indent=2))
    return str(output)
