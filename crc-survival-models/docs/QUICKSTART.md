# Quickstart / 快速开始

## 1. Clone the lightweight project / 克隆轻量项目

This public subproject contains only the CRC survival models and their inference code.
UNI and CellViT++ source/checkpoints are not redistributed here.

公开子项目只包含 CRC 生存模型及其推理代码，不再分发 UNI 和 CellViT++ 源码或权重：

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/crc-survival-models
python predict.py validate-install
```

Official upstream repositories / 官方原库：

- CellViT++: https://github.com/TIO-IKIM/CellViT-plus-plus
- UNI: https://github.com/mahmoodlab/UNI

Prepared-feature mode needs no third-party source checkout. For raw WSI mode, obtain
compatible files under the upstream terms and either pass the four path arguments shown
below or set `CRC_SURVIVAL_CELLVIT_ROOT`, `CRC_SURVIVAL_CELLVIT_CHECKPOINT`,
`CRC_SURVIVAL_CELLVIT_CLASSIFIERS`, and `CRC_SURVIVAL_UNI_CHECKPOINT`.

已有特征模式无需第三方源码。原始 WSI 模式需按官方许可自行获取兼容文件，并通过
下方四个参数或对应环境变量提供本地路径。项目不会自动下载或重新分发这些文件。

## 2. Create the environment / 创建环境

```bash
conda env create -f environment.yml
conda activate crc-survival
```

This command is documentation for the future downloader; no package was installed while preparing this staging folder. / 这是给未来下载者的安装说明；本次整理过程中没有擅自安装任何依赖。

## 3. Preflight a WSI / 预检切片

```bash
python predict_wsi.py slide.svs --case-id CRC_001 --model both \
  --output-dir results/CRC_001 --check-only \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

Alternatively, after setting the four environment variables, run
`python predict.py validate-install --raw-wsi`. / 也可先设置四个环境变量，再运行
`python predict.py validate-install --raw-wsi` 检查原始 WSI 依赖。

If the WSI lacks valid physical-resolution metadata, add `--mpp 0.25` only when that value is known from the scanner or acquisition record. Never guess MPP. / 若切片缺少物理分辨率元数据，仅在扫描记录明确时使用 `--mpp 0.25`；禁止猜测 MPP。

## 4. Run one or several slides / 运行单张或多张切片

```bash
python predict_wsi.py slide_1.svs slide_2.svs \
  --case-id CRC_001 --model both --device cuda:0 \
  --output-dir results/CRC_001 \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

Output / 输出：

- `survival_predictions.json`: both final continuous risks and frozen groups.
- `cell_topology/`: typed nuclei, cell graph, topology tokens and topology prediction.
- `patch_feature/`: patch cache, extraction audit and Patch-model prediction.

The command is resumable at completed intermediate files. Keep all slides for one patient under the same `case-id`. / 完成的中间文件可复用；同一患者的所有切片必须使用同一个 `case-id`。

## 5. Prepared-feature mode / 已有特征模式

Cell-topology token files:

```bash
python predict.py topology *.topology_tokens.h5 --case-id CRC_001 \
  --output results/CRC_001/topology.json
```

Patch caches:

```bash
python predict.py patch *.patch_cache.h5 --case-id CRC_001 \
  --device cuda:0 --output results/CRC_001/patch.json
```

## 6. Quality-control failures / 质控失败

The pipeline stops instead of fabricating a score when:

- slide metadata or physical scale is invalid;
- native resolution is too low for nuclear analysis;
- tissue detection is empty;
- too few valid tumour-front topology regions are available;
- an HDF5 feature contract is incomplete or inconsistent.

遇到以下情况流程会停止，而不会强行给出风险值：物理尺度无效、原始分辨率过低、组织检测为空、肿瘤前沿区域不足、或 HDF5 特征契约不完整。

## 7. Interpretation / 结果解释

Compare patients only after identical preprocessing and QC. Use the continuous score for modelling. The zero-threshold group is a frozen descriptive label, not a calibrated clinical decision boundary. Do not reuse Fig3 outcome-adaptive cutpoints on new patients.

患者间比较必须使用完全一致的预处理和质控。统计建模优先使用连续分数；零阈值分组只是冻结的描述性标签，不是临床决策界值。不得把 Fig3 的结局自适应切点用于新患者。
