# CRC Survival Models

两套结直肠癌（CRC）病理生存风险模型：输入全切片图像（WSI）或已提取的特征文件，
输出连续的总体生存相对风险。`Fig3` 生存曲线由这两套模型生成。

| 模型 | 输入 | 主要思路 |
|---|---|---|
| Cell Topology | 带细胞类型的空间邻接图 | 比较全组织和肿瘤前沿的细胞拓扑，不读取图像颜色、原始形态、临床变量或微生物数据 |
| Patch Feature | 病理 Patch 表征及空间邻域 | 联合局部表征梯度、肿瘤核心、侵袭前沿和瘤周结构，由五折模型给出患者风险 |

两套模型都输出连续的总体生存相对风险。`risk >= 0` 记为高风险，`risk < 0`
记为低风险。这个分组只是固定阈值下的描述，不是生存概率、剩余寿命或治疗建议。

## 快速开始

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/crc-survival-models
conda env create -f environment.yml
conda activate crc-survival
python predict.py validate-install
```

已有特征文件时，直接打分：

```bash
python predict.py topology slide_a.topology_tokens.h5 slide_b.topology_tokens.h5 \
  --case-id CRC_001 --output results/topology.json
```

`validate-install` 只检查本项目自带的生存模型。原始 WSI 流程还需要用户在本地准备
CellViT++ 和 UNI 相关文件，见[外部依赖说明](docs/THIRD_PARTY_NOTICES.md)。本仓库不复制
这两个项目的源码或权重。

## 用法

### 已有特征文件

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

两类输入文件由本仓库的 WSI 流程生成，字段约定见
[Cell Topology 模型卡](models/cell_topology/MODEL_CARD.md) 与
[Patch Feature 模型卡](models/patch_feature/MODEL_CARD.md)。同一患者有多张切片时，
应在同一条命令中输入，不能把每张切片当成独立患者。

### 原始 WSI

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

完整的 WSI 流程与失败停止条件见 [WSI 流程说明](docs/QUICKSTART.md)。

### Python 调用

两个模型都提供 `predict()`，返回值与 `predict.py` 写出的 JSON 一致。在仓库根目录运行：

```python
from models.cell_topology.model import predict

result = predict(
    ["slide_a.topology_tokens.h5", "slide_b.topology_tokens.h5"],
    case_id="CRC_001",
)
```

```python
from models.patch_feature.model import predict

result = predict(
    ["slide_a.patch_cache.h5", "slide_b.patch_cache.h5"],
    case_id="CRC_001",
    device="cuda:0",  # 省略时自动选择 GPU 或 CPU
)
```

两个函数都接受可选的权重目录参数（细胞拓扑为 `model_dir`，Patch 为 `weights_dir`），
默认使用各自 `weights/` 目录中的随附权重。

## 输出说明

`predict.py topology` 输出：

| 字段 | 说明 |
|---|---|
| `model` | `"cell_topology"` |
| `case_id`、`n_slides` | 病例标识；参与打分的切片数 |
| `risk` | 连续相对风险分数（全组织与肿瘤前沿两个分支的等权平均） |
| `risk_group`、`threshold` | `risk >= 0` 为 `"high"`，否则 `"low"`；阈值固定为 0 |
| `branches` | `whole_tissue` 与 `tumour_front` 分支各自的 `risk` 与 `n_regions` |

`predict.py patch` 输出：

| 字段 | 说明 |
|---|---|
| `model` | `"patch_feature"` |
| `case_id`、`n_slides` | 病例标识；参与打分的切片数 |
| `risk` | 连续相对风险分数（五折风险的中位数） |
| `risk_group`、`threshold` | 同上 |
| `fold_mean`、`fold_sd` | 五折风险的均值与标准差 |
| `fold_predictions` | 五折各自的 `raw_risk` 与标准化 `risk` |

`predict_wsi.py` 在 `--output-dir` 下写出：

```text
results/CRC_001/
├── survival_predictions.json    # case_id、slides 与两套模型的完整结果
├── cell_topology/               # 切片清单、细胞检测状态、细胞核与图、拓扑 token
└── patch_feature/               # 切片清单、patch cache 与提取记录
```

## 测试与验证

- 结构检查：`python -m pytest tests/ -q`
- 安装校验：`python predict.py validate-install`，加 `--raw-wsi` 可额外检查本地
  CellViT++ 与 UNI 路径
- 验证记录：`docs/validation/` 保留两次真实运行的完整数值，可用来核对运行环境：
  - 304 MB 结直肠癌 SVS 端到端：细胞拓扑 risk `-0.3449`（low）、Patch risk `0.2766`
    （high）。两套模型使用不同输入，同一病例可以落入不同分组。
  - Patch 模型 CPU 复跑与 GPU 记录值的连续分数差约 0.0006，分组一致；这一量级的
    差异来自浮点执行环境。

## Fig3 结果

细胞拓扑模型的 Fig3 结果来自嵌套五折 OOF 预测：431 名患者、92 个事件，
C-index 0.6454，log-rank P=`8.75e-05`，每 IQR 风险增量 HR 1.810
（95% CI 1.329 至 2.464）。

Patch 模型在四个队列上的连续风险随机效应 Meta 分析为 HR 1.743
（95% CI 1.176 至 2.583），P=`0.0206`。各队列统计量见
[Fig3 模型溯源](docs/FIG3_PROVENANCE.md)。

## 项目结构

```text
crc-survival-models/
├── predict.py                   已有特征文件的推理入口（topology / patch / validate-install）
├── predict_wsi.py               原始 WSI 端到端入口（含 --check-only 预检）
├── environment.yml              conda 环境（Python 3.10、PyTorch 2.2.1、OpenSlide 等）
├── models/
│   ├── cell_topology/           模型卡、推理代码、权重、WSI 特征流程、参考统计量
│   └── patch_feature/           模型卡、推理代码、五折权重、Patch cache 提取器、参考统计量
├── docs/
│   ├── QUICKSTART.md            WSI 完整流程
│   ├── FIG3_PROVENANCE.md       Fig3 与参考统计量的对应关系
│   ├── THIRD_PARTY_NOTICES.md   外部依赖来源与许可
│   └── validation/              验证运行记录
├── tests/                       结构检查（pytest）
├── SKILL.md                     供 Agent 调用本项目时的操作约定
└── LICENSE                      CC BY-NC-SA 4.0
```

## 限制与适用范围

模型仅供科研使用，不用于诊断、治疗选择或个体生存期估计。详细限制见两份模型卡：
[Cell Topology](models/cell_topology/MODEL_CARD.md) 和
[Patch Feature](models/patch_feature/MODEL_CARD.md)。

## 许可与引用

本子项目采用 **CC BY-NC-SA 4.0**（署名—非商业性使用—相同方式共享）许可，完整文本见
[LICENSE](LICENSE)；许可只覆盖仓库原创代码、文档和随附模型，不改变用户另行取得的
UNI 或 CellViT++ 材料的许可证。使用请注明来源并引用本仓库；论文信息发表后在此补充。

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
probability or a clinical threshold.

Prepared topology tokens and patch caches can be scored with `predict.py`. Raw WSI
inference uses `predict_wsi.py` and requires locally supplied compatible CellViT++ and
UNI-related files. This repository does not redistribute those upstream projects or
their checkpoints. Both models also expose `predict()` for use from Python.

For the Cell Topology model, the nested five-fold OOF result was n=431 with 92 events,
C-index 0.6454, log-rank P=`8.75e-05`, and HR 1.810 per IQR of risk
(95% CI 1.329 to 2.464). For the Patch Feature model, the continuous-risk random-effects
estimate across four cohorts was HR 1.743 (95% CI 1.176 to 2.583), P=`0.0206`.

Structure checks run with `python -m pytest tests/ -q`; validation records from two full
runs are kept under `docs/validation/`. This code and the included models are for
research use only, and are released under CC BY-NC-SA 4.0 (see `LICENSE`); please cite
this repository when using them.

</details>
