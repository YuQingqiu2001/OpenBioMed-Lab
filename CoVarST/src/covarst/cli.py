"""CoVarST 公开入口；冻结研究脚本保留原始参数与帮助。"""
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

COMMAND_HELP={'list-models': '列出癌种、参考规模、教师数量及模型参数量', 'verify-assets': '逐文件校验模型和参考资产的 SHA256', 'infer-spots': 'H&E 特征推理为 spot 全基因表达概率', 'extract-wsi': '从真实 H&E 切片提取 ParamNet/Virchow2 特征', 'prepare-coarse': '准备同切片实测粗 ST 与真实细胞核输入', 'fit-coarse': '拟合粗 ST 到细胞状态；详细参数见 fit-coarse -- --help', 'materialize-cells': '分块生成全基因归一化细胞表达概率 P_post', 'allocate-cell-mass': '生成守恒分数计数 X_mass 与未分配质量', 'prepare-reference-counts': '整理原始计数及癌种基因面板', 'deconvolve-local': '逐切片原版 BayesTME 分解，需单独兼容环境', 'compress-programs': '早期 RCCD 辅助入口，默认值不是最终方案', 'fit-fixed-reference': '按锁定划分拟合并选择最终共识参考', 'refine-reference': '联合细化参考、spot theta 与批次适配器', 'train-he': '在固定参考下训练原架构图像映射器', 'prepare-native16': '准备 native16 核普查与粗 RNA 谱', 'prepare-virtual55': '按真实捕获几何准备 virtual55 降采样', 'fit-native16': '用严格内存受限投影器拟合 native16', 'fit-virtual55': '用连通分量测量零空间投影器拟合 virtual55'}
ARGUMENT_HELP={'--assets': 'CoVarST 资产根目录，需包含 models/ 与 references/', '--cancer': '准确癌种标识，例如 colorectal_cancer', '--features': '图像特征 NPZ，含 features、coords、barcode', '--output': '新输出路径，不覆盖既有结果', '--device': 'cpu 或 cuda；默认值见使用指南', '--batch': '精确登记的 SOURCE|TECHNOLOGY，省略时使用中性参考', '--gene-block': '解码基因块大小，默认 512', '--wsi': '真实 H&E 切片文件', '--paramnet-root': '用户提供的 ParamNet 源码与权重目录', '--virchow2-root': '用户提供的 Virchow2 权重目录', '--mpp': '有实测依据的微米/像素标度，省略时读取元数据', '--batch-size': '图块批量，默认 16', '--input': '同切片实测粗计数与真实核的规范 HDF5', '--sample': '切片标识', '--seed': '随机种子，默认 20261010', '--factors': '已完成的细胞因子 HDF5', '--cell-block': '细胞块大小，默认 256', '--prepared': '准备完成的粗 ST HDF5', '--help': '显示帮助并退出'}

def main(argv=None):
 argv=list(sys.argv[1:] if argv is None else argv)
 if argv and argv[0] in RECIPES:
  name=argv.pop(0);engine,script=RECIPES[name];activate('he','reference','reference/upstream_pipeline',engine)
  if argv[:1]==['--']:argv.pop(0)
  if name=='fit-virtual55':argv=['--engine','v55']+argv
  sys.argv=[str(ENGINE/engine/script)]+argv;runpy.run_path(sys.argv[0],run_name='__main__');return
 p=argparse.ArgumentParser(description='CoVarST：每癌种一个 H&E 模型，同切片粗 ST 重建细胞表达');sub=p.add_subparsers(dest='command',required=True)
 for command in ['list-models','verify-assets']:
  q=sub.add_parser(command,help=COMMAND_HELP[command]);q.add_argument('--assets',type=Path)
 q=sub.add_parser('infer-spots',help=COMMAND_HELP['infer-spots']);q.add_argument('--cancer',required=True);q.add_argument('--features',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--assets',type=Path);q.add_argument('--device',default='cpu');q.add_argument('--batch');q.add_argument('--gene-block',type=int,default=512)
 q=sub.add_parser('extract-wsi',help=COMMAND_HELP['extract-wsi']);q.add_argument('--wsi',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--paramnet-root',type=Path,required=True);q.add_argument('--virchow2-root',type=Path,required=True);q.add_argument('--mpp',type=float);q.add_argument('--device',default='cuda');q.add_argument('--batch-size',type=int,default=16)
 q=sub.add_parser('prepare-coarse',help=COMMAND_HELP['prepare-coarse']);q.add_argument('--input',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--sample',required=True);q.add_argument('--seed',type=int,default=20261010)
 q=sub.add_parser('fit-coarse',help=COMMAND_HELP['fit-coarse']);q.add_argument('arguments',nargs=argparse.REMAINDER)
 q=sub.add_parser('materialize-cells',help=COMMAND_HELP['materialize-cells']);q.add_argument('--factors',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--cell-block',type=int,default=256);q.add_argument('--gene-block',type=int,default=512)
 q=sub.add_parser('allocate-cell-mass',help=COMMAND_HELP['allocate-cell-mass']);q.add_argument('--prepared',type=Path,required=True);q.add_argument('--factors',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 for name in RECIPES:sub.add_parser(name,help=COMMAND_HELP[name])
 p._positionals.title='命令';p._optionals.title='可选参数'
 for parser in [p]+list(sub.choices.values()):
  parser._positionals.title='命令或位置参数';parser._optionals.title='可选参数'
  for action in parser._actions:
   for option in action.option_strings:
    if option in ARGUMENT_HELP:action.help=ARGUMENT_HELP[option];break
   if action.dest=='arguments':action.help='通过 -- 转交冻结拟合器的参数'
 for name,parser in sub.choices.items():parser.description=COMMAND_HELP[name]
 a=p.parse_args(argv)
 if a.command in ('list-models','verify-assets'):
  root=asset_root(a.assets);file=root/'models/manifest.json'
  if not file.exists():raise FileNotFoundError('未找到最终模型发布清单')
  manifest=json.loads(file.read_text());
  if a.command=='verify-assets':
   for item in manifest['files']:
    if sha256(root/item['file'])!=item['sha256']:raise ValueError('文件哈希不一致：'+item['file'])
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
