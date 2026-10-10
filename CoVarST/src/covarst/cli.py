"""CoVarST public entry points; frozen recipes are exposed with their own help."""
import argparse, json,runpy,sys
from pathlib import Path
from .runtime import activate,asset_root,ENGINE,sha256

RECIPES={
 'prepare-reference-counts':('reference','prepare_hest_tissue_top10k.py'),
 'deconvolve-local':('reference','upstream_pipeline/run_pure_upstream_bayestme_per_slide.py'),
 'compress-programs':('reference','upstream_pipeline/fit_bayestme_top10k_rccd_consensus.py'),
 'fit-fixed-reference':('reference','upstream_pipeline/fit_bayestme_top10k_direct_consensus_locked_split.py'),
 'refine-reference':('reference','upstream_pipeline/refine_direct_consensus_spot_nmf.py'),
 'train-he':('he','train_virchow2_pan_atlas_consensus_mapper.py'),
 'prepare-native16':('cells','prepare_generic_hd_inputs.py'),
 'prepare-virtual55':('cells','prepare_generic_v55_inputs.py'),
 'fit-native16':('cells','fit_bayesgraph_qr_memory_v3.py'),
 'fit-virtual55':('cells','fit_bayesgraph_component_projection.py')}

def main(argv=None):
 argv=list(sys.argv[1:] if argv is None else argv)
 if argv and argv[0] in RECIPES:
  name=argv.pop(0);engine,script=RECIPES[name];activate('he','reference','reference/upstream_pipeline',engine)
  if argv[:1]==['--']:argv.pop(0)
  if name=='fit-virtual55':argv=['--engine','v55']+argv
  sys.argv=[str(ENGINE/engine/script)]+argv;runpy.run_path(sys.argv[0],run_name='__main__');return
 p=argparse.ArgumentParser(description='CoVarST: one deployable H&E checkpoint per cancer');sub=p.add_subparsers(dest='command',required=True)
 for command in ['list-models','verify-assets']:
  q=sub.add_parser(command);q.add_argument('--assets',type=Path)
 q=sub.add_parser('infer-spots');q.add_argument('--cancer',required=True);q.add_argument('--features',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--assets',type=Path);q.add_argument('--device',default='cpu');q.add_argument('--batch');q.add_argument('--gene-block',type=int,default=512)
 q=sub.add_parser('extract-wsi');q.add_argument('--wsi',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--paramnet-root',type=Path,required=True);q.add_argument('--virchow2-root',type=Path,required=True);q.add_argument('--mpp',type=float);q.add_argument('--device',default='cuda');q.add_argument('--batch-size',type=int,default=16)
 q=sub.add_parser('prepare-coarse');q.add_argument('--input',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--sample',required=True);q.add_argument('--seed',type=int,default=20261010)
 q=sub.add_parser('fit-coarse');q.add_argument('arguments',nargs=argparse.REMAINDER)
 q=sub.add_parser('materialize-cells');q.add_argument('--factors',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--cell-block',type=int,default=256);q.add_argument('--gene-block',type=int,default=512)
 q=sub.add_parser('allocate-cell-mass');q.add_argument('--prepared',type=Path,required=True);q.add_argument('--factors',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 for name in RECIPES:sub.add_parser(name,help='Frozen upstream recipe; invoke COMMAND --help for arguments')
 a=p.parse_args(argv)
 if a.command in ('list-models','verify-assets'):
  root=asset_root(a.assets);file=root/'models/manifest.json'
  if not file.exists():raise FileNotFoundError('Final checkpoint manifest has not been built yet')
  manifest=json.loads(file.read_text());
  if a.command=='verify-assets':
   for item in manifest['files']:
    if sha256(root/item['file'])!=item['sha256']:raise ValueError('Hash mismatch: '+item['file'])
   print(json.dumps({'verified':True,'cancers':len(manifest['models']),'files':len(manifest['files'])}));return
  print(json.dumps(manifest['models'],ensure_ascii=False,indent=2));return
 if a.command=='infer-spots':
  from .he import infer_spots
  print(json.dumps(infer_spots(a.cancer,a.features,a.output,a.assets,a.device,a.batch,a.gene_block)));return
 if a.command=='extract-wsi':
  from .he import extract_wsi
  print(extract_wsi(a.wsi,a.output,a.paramnet_root,a.virchow2_root,a.device,a.mpp,a.batch_size));return
 from . import cells
 if a.command=='prepare-coarse':print(json.dumps(cells.prepare_coarse(a.input,a.output,a.sample,a.seed)))
 elif a.command=='fit-coarse':cells.fit_coarse(a.arguments)
 elif a.command=='materialize-cells':print(json.dumps(cells.materialize(a.factors,a.output,a.cell_block,a.gene_block)))
 elif a.command=='allocate-cell-mass':print(json.dumps(cells.allocate_mass(a.prepared,a.factors,a.output)))
