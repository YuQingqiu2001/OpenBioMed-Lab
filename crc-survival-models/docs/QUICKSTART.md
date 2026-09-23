# 快速开始

## 1. 建立环境

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/crc-survival-models
conda env create -f environment.yml
conda activate crc-survival
python predict.py validate-install
```

这一步只检查仓库自带的两套生存模型。CellViT++ 和 UNI 相关文件不在仓库中。

## 2. 已有特征文件时直接计算风险

细胞拓扑 token：

```bash
python predict.py topology *.topology_tokens.h5 \
  --case-id CRC_001 --output results/CRC_001/topology.json
```

Patch cache：

```bash
python predict.py patch *.patch_cache.h5 \
  --case-id CRC_001 --device cuda:0 \
  --output results/CRC_001/patch.json
```

同一患者的所有切片应共用一个 `case-id`，并在一次调用中输入。

## 3. 配置原始 WSI 流程

从官方项目取得所需文件：

- [CellViT++](https://github.com/TIO-IKIM/CellViT-plus-plus)
- [UNI](https://github.com/mahmoodlab/UNI)

调用时可以直接给出四个路径：

```bash
python predict_wsi.py slide.svs \
  --case-id CRC_001 --model both \
  --output-dir results/CRC_001 --check-only \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

也可以设置以下环境变量，再运行 `python predict.py validate-install --raw-wsi`：

```text
CRC_SURVIVAL_CELLVIT_ROOT
CRC_SURVIVAL_CELLVIT_CHECKPOINT
CRC_SURVIVAL_CELLVIT_CLASSIFIERS
CRC_SURVIVAL_UNI_CHECKPOINT
```

官方仓库不一定包含本流程所需的项目特定分类头。文件来源和兼容性要求见
[外部依赖说明](THIRD_PARTY_NOTICES.md)。

## 4. 运行 WSI

```bash
python predict_wsi.py slide_1.svs slide_2.svs \
  --case-id CRC_001 --model both --device cuda:0 \
  --output-dir results/CRC_001 \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

如果切片缺少物理分辨率信息，只有在扫描记录明确时才使用 `--mpp`。不要根据文件名或
经验猜测 MPP。

程序会复用已经完成的中间文件。主要输出为：

- `survival_predictions.json`：两套模型的连续风险和固定阈值分组
- `cell_topology/`：细胞核、细胞图、拓扑 token 和拓扑风险
- `patch_feature/`：Patch cache、提取记录和 Patch 风险

## 5. 什么时候会停止

遇到以下情况时，流程不会给出风险值：

- 切片缺少可靠的物理尺度
- 原始分辨率不足以进行细胞核分析
- 没有检测到有效组织
- 肿瘤前沿区域数量不足
- HDF5 文件不符合模型输入约定

不同患者之间的风险只能在相同预处理和质控条件下比较。

<details>
<summary>English notes</summary>

Create the environment from `environment.yml`, then run
`python predict.py validate-install`. Prepared topology tokens and patch caches can be
scored directly with `predict.py`.

Raw WSI inference needs local compatible CellViT++ source, a segmentation checkpoint,
seven phenotype classifier heads, and a compatible UNI-derived coarse classifier.
Provide them with the four command-line options above or the corresponding environment
variables. Run with `--check-only` before starting a full slide.

All slides from one patient belong under one `case-id`. Supply `--mpp` only when the
scanner or acquisition record provides the value. The program stops on invalid scale,
insufficient resolution or tissue, too few tumour-front regions, and incompatible HDF5
inputs.

</details>
