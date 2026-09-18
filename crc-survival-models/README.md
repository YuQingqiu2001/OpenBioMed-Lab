# CRC Survival Models

这里保存的是 `Fig3` 生存曲线实际使用的两套结直肠癌病理模型。

| 模型 | 输入 | 主要思路 |
|---|---|---|
| Cell Topology | 带细胞类型的空间邻接图 | 比较全组织和肿瘤前沿的细胞拓扑，不读取图像颜色、原始形态、临床变量或微生物数据 |
| Patch Feature | 病理 Patch 表征及空间邻域 | 联合局部表征梯度、肿瘤核心、侵袭前沿和瘤周结构，由五折模型给出患者风险 |

两套模型都输出连续的总体生存相对风险。`risk >= 0` 记为高风险，`risk < 0`
记为低风险。这个分组只是固定阈值下的描述，不是生存概率、剩余寿命或治疗建议。

## 安装

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/crc-survival-models
conda env create -f environment.yml
conda activate crc-survival
python predict.py validate-install
```

`validate-install`只检查本项目自带的生存模型。原始 WSI 流程还需要用户在本地准备
CellViT++ 和 UNI 相关文件，见[外部依赖说明](docs/THIRD_PARTY_NOTICES.md)。本仓库不复制
这两个项目的源码或权重。

## 已有特征文件

细胞拓扑 token：

```bash
python predict.py topology slide_a.topology_tokens.h5 slide_b.topology_tokens.h5 \
  --case-id CRC_001 --output results/topology.json
```

Patch cache：

```bash
python predict.py patch slide_a.patch_cache.h5 slide_b.patch_cache.h5 \
  --case-id CRC_001 --device cuda:0 --output results/patch.json
```

同一患者有多张切片时，应在同一条命令中输入，不能把每张切片当成独立患者。

## 从 WSI 开始

```bash
python predict_wsi.py patient_slide.svs \
  --case-id CRC_001 --model both --output-dir results/CRC_001 \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

支持 `SVS`、`MRXS`、`NDPI`、`TIF/TIFF`、`SCN` 和 `BIF`。Patch 分支还支持项目
自定义的 `IBL` 读取器。正式运行前先加上 `--check-only`，检查切片、第三方文件和模型路径。

细胞拓扑分支建议使用约 24 GB 显存的 CUDA GPU 和至少 32 GB 内存。原始切片的分辨率
不能低于 0.75 µm/pixel。Windows 用户宜在 WSL2 或 Linux 中运行 WSI 流程；已有特征
文件的推理没有这个限制。

## Fig3 与新样本的阈值

Fig3 中的 Patch 曲线采用了各队列基于结局选择的切点，因此分离程度会受到选择偏倚影响。
这些切点保留在 `reference_results` 中用于复现原图，不应用于新患者。新样本应报告连续风险，
如需描述性分组，则使用固定阈值 0。

细胞拓扑模型的 Fig3 结果来自嵌套五折 OOF 预测：431 名患者、92 个事件，
C-index 0.6454，log-rank P=`8.75e-05`，每 IQR 风险增量 HR 1.810
（95% CI 1.329 至 2.464）。

Patch 模型在四个队列上的连续风险随机效应 Meta 分析为 HR 1.743
（95% CI 1.176 至 2.583），P=`0.0206`。图中切点及各队列结果见
[Fig3 模型溯源](docs/FIG3_PROVENANCE.md)。

## 仓库内容

- `models/cell_topology/`：细胞拓扑模型、WSI 特征流程和参考统计量
- `models/patch_feature/`：五折 Patch 模型、Patch cache 提取器和参考统计量
- `predict.py`：已有特征文件的推理入口
- `predict_wsi.py`：原始 WSI 推理入口
- `docs/`：使用说明、模型与统计量对应关系、外部依赖说明
- `SKILL.md`：供 Agent 调用本项目时使用的操作约定

模型仅供科研使用，不用于诊断、治疗选择或个体生存期估计。详细限制见两份模型卡：
[Cell Topology](models/cell_topology/MODEL_CARD.md) 和
[Patch Feature](models/patch_feature/MODEL_CARD.md)。

<details>
<summary>English summary</summary>

This subproject contains the two colorectal cancer pathology models used for the
survival curves in `Fig3`.

The Cell Topology model uses typed-cell neighbourhood graphs from whole tissue and the
tumour front. It does not use image colour, raw morphology, clinical variables, or
microbial measurements. The Patch Feature model combines histology patch embeddings,
local spatial gradients, and soft tumour-core, invasive-front, and peritumour regions.

Both models return a continuous relative-risk score for overall survival. The fixed
descriptive split is `risk >= 0` versus `risk < 0`; it is not a calibrated survival
probability or a clinical threshold. The outcome-adaptive cutpoints shown in some Fig3
panels are provided only for figure reproduction and must not be reused for new patients.

Prepared topology tokens and patch caches can be scored with `predict.py`. Raw WSI
inference uses `predict_wsi.py` and requires locally supplied compatible CellViT++ and
UNI-related files. This repository does not redistribute those upstream projects or
their checkpoints.

For the Cell Topology model, the nested five-fold OOF result was n=431 with 92 events,
C-index 0.6454, log-rank P=`8.75e-05`, and HR 1.810 per IQR of risk
(95% CI 1.329 to 2.464). For the Patch Feature model, the continuous-risk random-effects
estimate across four cohorts was HR 1.743 (95% CI 1.176 to 2.583), P=`0.0206`.

This code and the included models are for research use only.

</details>
